import numpy as np
import matplotlib.pyplot as plt
import dynamics.vibration as dy

acceptance_level = np.array([[  20, 0.0053],
                             [ 150, 0.04],
                             [ 600, 0.04],
                             [2000, 0.0036]])



freq_q = np.linspace(20, 2000, 1000)
sdof_response = dy.sdof_psd_response(base_spectrum=acceptance_level, f_n=np.array([100, 200, 300]), Q=10, f_query=freq_q)
vrs, f_vrs = dy.vrs(acceptance_level)
vrs_nsigma, f_vrs_nsigma = dy.vrs_shock_equivalent(acceptance_level, 60)
vrs_miles, f_vrs_miles = dy.vrs_miles(spectrum=acceptance_level, Q=10, duration=60)


fig, ax = plt.subplots()
plt.loglog((acceptance_level[:,0]), (acceptance_level[:,1]))
plt.loglog(freq_q, sdof_response.T, "--")
ax.set_ylim(0.001, 10)
ax.grid()
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='gray') 
plt.xlabel("Frequency (Hz)")
plt.ylabel(r'Acceleration PSD ($g^2$/Hz)')
plt.title('SDOF Response to MIL-STD-1540C ATP Base Input\nSDOF Systems $f_n$ = 100, 200, 300 Hz')
plt.savefig("sdof_response_python.png", dpi=300)


fig, ax2 = plt.subplots()
plt.loglog(f_vrs, vrs)
ax2.set_ylim(1, 20)
ax2.grid()
ax2.grid(which='minor', linestyle=':', linewidth='0.5', color='gray') 
plt.xlabel("Natural Frequency (Hz)")
plt.ylabel("Acceleration (gRMS)")
plt.title('Vibration Response Spectrum\nSDOF systems Q=10   Base Input: MIL-STD-1540C ATP')
plt.savefig("vrs_python.png", dpi=300)


fig, ax3 = plt.subplots()
plt.loglog(f_vrs_miles, vrs_miles)
plt.loglog(f_vrs_nsigma, vrs_nsigma)
ax3.grid()
ax3.grid(which='minor', linestyle=':', linewidth='0.5', color='gray') 
plt.xlabel("Natural Frequency (Hz)")
plt.ylabel(r"Acceleration (g)")
plt.legend(["VRS, Miles Eq. Approximation", "VRS, shock equivalence"])
plt.title('Vibration Response Spectrum, Q=10, SRS Equivalence\nMIL-STD-1540C ATP 60 sec duration')
plt.savefig("vrs_nsigma_python.png", dpi=300)