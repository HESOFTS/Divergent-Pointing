import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
# import tempfile
from copy import deepcopy


import pandas as pd
from astropy.coordinates.angle_utilities import angular_separation

from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time

from ctapipe.io import event_source
from ctapipe.io import HDF5TableWriter
from ctapipe.calib import CameraCalibrator
from ctapipe.reco import HillasReconstructor
from ctapipe.image import hillas_parameters, leakage
from ctapipe.image.cleaning import tailcuts_clean
# from ctapipe.visualization import ArrayDisplay
# from ctapipe.image.timing_parameters import timing_parameters
import tables
import iof
import sys
import os

# Import a simtelarray file of MC simulated data
# input_file = "/Users/alicedonini/PycharmProjects/ctapipe/gamma_20deg_180deg_run875___cta-prod3-demo-2147m-LaPalma-baseline.simtel.gz"
input_file = sys.argv[1]
output_filename = ('dl1_' + os.path.basename(input_file).split('.')[0] + '.h5')




events = event_source(input_file, back_seekable=True, max_events=None)

calibrator = CameraCalibrator()

location = EarthLocation.of_site('Roque de los Muchachos')
obstime = Time('2018-11-01T02:00')
horizon_frame = AltAz(location=location, obstime=obstime)

cleaning_level = {
    'LSTCam': (4, 8, 2),
    'NectarCam': (5, 10, 2),
    'FlashCam': (4, 8, 2),
    'std': (4, 8, 2)
}



filters = tables.Filters(
    complevel=5,    # enable compression, with level 0=disabled, 9=max
    complib='blosc:zstd',   #  compression using blosc
    fletcher32=True,    # attach a checksum to each chunk for error correction
    bitshuffle=False,   # for BLOSC, shuffle bits for better compression
)




### Write extra information to the DL1 file
event = next(iter(events))
metadata = iof.global_metadata(events)
iof.write_metadata(metadata, output_filename)
iof.write_array_info(event, output_filename)
iof.write_mcheader(event.mcheader, output_filename, obs_id=event.r0.obs_id, filters=filters, metadata=metadata)


extra_im_c = iof.ExtraImageInfo()
event_info_c = iof.EventInfo()

with HDF5TableWriter(filename=output_filename, group_name='events', mode='a') as writer:

    for ii, event in enumerate(events):
        if ii%1 == 0:
            print(ii)
        calibrator(event)
        if len(event.r0.tels_with_data) < 2:
            continue

        reco = HillasReconstructor()

        # mapping of telescope_id to parameters for stereo reconstruction
        hillas_dict = {}
        #time_gradients = {}

        telescope_pointings = {}

        for telescope_id, dl1 in event.dl1.tel.items():

            tel = event.dl1.tel[telescope_id]
            tel.prefix = ''  # don't really need one for HDF5 writer
            extra_im_c.tel_id = telescope_id
            extra_im_c.event_id = event.r0.event_id
            extra_im_c.obs_id = event.r0.obs_id

            camera = event.inst.subarray.tels[telescope_id].camera
            image = dl1.image
            #pulse_time = dl1.pulse_time

            if camera.cam_id in cleaning_level:
                boundary, picture, min_neighbors = cleaning_level[camera.cam_id]
            else:
                boundary, picture, min_neighbors = cleaning_level['std']

            # create a clean mask of pixels above the threshold
            clean = tailcuts_clean(
                camera,
                image,
                boundary_thresh=boundary,
                picture_thresh=picture,
                min_number_picture_neighbors=min_neighbors
            )

            # require more than five pixels after cleaning in each telescope
            if clean.sum() < 5:
                continue

            hillas_c = hillas_parameters(camera[clean], image[clean])

            leakage_c = leakage(camera, image, clean)
            # n_islands, island_ids = number_of_islands(camera, clean)

            # remove if after cleaning we have more than 1 island
            # if n_islands > 1:
            #     continue

            # remove events with high leakage
            if leakage_c.leakage2_intensity > 0.2:
                continue

            hillas_dict[telescope_id] = hillas_c
            #timing_c = timing_parameters(camera[clean], image[clean], pulse_time[clean], hillas_c)
            # ssts have no timing in prod3b, so we'll use the skewness
            #time_gradients[telescope_id] = timing_c.slope.value if camera.cam_id != 'ASTRICam' else hillas_c.skewness

            # this makes sure, that we get an arrow in the array plow for each telescope
            # might have the wrong direction though
            #if abs(time_gradients[telescope_id]) < 0.2:
            #    time_gradients[telescope_id] = 1.0

            # Divergent part
            telescope_pointings[telescope_id] = SkyCoord(alt=event.mc.tel[telescope_id].altitude_raw * u.rad,
                                                         az=event.mc.tel[telescope_id].azimuth_raw * u.rad,
                                                         frame=horizon_frame)

            # remove the first part of the tel_name which is the type 'LST', 'MST' or 'SST'
            tel_name = str(event.inst.subarray.tel[telescope_id])[4:]
            writer.write(table_name='dl1', containers=[hillas_c, leakage_c, extra_im_c])
            
            # Camera display

            # from ctapipe.visualization import CameraDisplay
            #
            # display = CameraDisplay(camera)
            # cleaned = dl1.image.copy()
            # cleaned[~clean] = 0.0
            #
            # display.image = cleaned
            # display.add_colorbar()
            #
            # display.overlay_moments(hillas_c, color='xkcd:red')
            # plt.figure()
            # plt.show()

        # array pointing needed for the creation of the TiltedFrame to perform the impact point reconstruction
        array_pointing = SkyCoord(az=event.mc.az, alt=event.mc.alt, frame=horizon_frame)

        event_info_c.event_id = event.r0.event_id
        event_info_c.obs_id = event.r0.obs_id

        if len(hillas_dict) > 1:
            stereo = reco.predict(
                hillas_dict, event.inst, array_pointing, telescope_pointings
            )

            writer.write(table_name='reconstructed',
                         containers=[stereo, event_info_c])
            writer.write(table_name='true', containers=[event.mc, event_info_c])

        # save event for plotting later
        # if event.count == 3:
        #     plotting_event = deepcopy(event)
        #     plotting_hillas = hillas_dict
        #     #plotting_timing = time_gradients
        #     plotting_stereo = stereo


# thrown the simtel energy histogram needs to loop over all the simtel file so we don't dump in case of tests
if events.max_events is None:
    iof.write_simtel_energy_histogram(events, output_filename, obs_id=event.dl0.obs_id, metadata=metadata)



# angle_offset = plotting_event.mcheader.run_array_direction[0]

# Array display
# disp = ArrayDisplay(plotting_event.inst.subarray)
# disp.set_vector_hillas(
#     plotting_hillas,
#     #time_gradient=plotting_timing,
#     angle_offset=angle_offset,
#     length=500
# )
#
# plt.scatter(
#     plotting_event.mc.core_x, plotting_event.mc.core_y,
#     s=200, c='k', marker='x', label='True Impact',
# )
# plt.scatter(
#     plotting_stereo.core_x, plotting_stereo.core_y,
#     s=200, c='r', marker='x', label='Estimated Impact',
# )
#
# plt.legend()
# plt.xlim(-400, 400)
# plt.ylim(-400, 400)





# theta-square plot
df_rec = pd.read_hdf(output_filename, key='events/reconstructed')
df_true = pd.read_hdf(output_filename, key='events/true')

# get angular offset between reconstructed shower direction and MC
# generated shower direction
theta = angular_separation(
    df_rec.az.values * u.deg, df_rec.alt.values * u.deg,
    df_true.az.values * u.deg, df_true.alt.values * u.deg,
)

plt.hist(theta.to(u.deg).value**2, bins=25, range=[0, 0.3])
plt.xlabel(r'$\theta2 / deg2$')
plt.ylabel("# of events")
plt.show()