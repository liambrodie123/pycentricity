#!/home/isobel.romero-shaw/anaconda3/bin/python3.6
import bilby as bb
import python_code.waveform as wf
import python_code.overlap as ovlp
import numpy as np


np.random.seed(54321)
outdir = '/home/isobel.romero-shaw/public_html/PYCENTRICITY/pycentricity/injection_recovery/pop_redundant_parameters'
label = 'injection'
# injection parameters
injection_parameters = dict(
    mass_1=35.0,
    mass_2=30.0,
    eccentricity=0.0,
    luminosity_distance=440.0,
    theta_jn=0.4,
    psi=0.1,
    phase=1.2,
    geocent_time=0.0,
    ra=45,
    dec=5.73,
    chi_1=0.0,
    chi_2=0.0,
)
# Frequency settings
minimum_frequency = 20
maximum_frequency = 1024
sampling_frequency = 4096
duration = 8
deltaT = 0.2
# Comparison waveform
comparison_waveform_generator = wf.get_IMRPhenomD_comparison_waveform_generator(
    minimum_frequency, sampling_frequency, duration
)
comparison_waveform_frequency_domain = comparison_waveform_generator.frequency_domain_strain(injection_parameters)
# Interferometers
interferometers = bb.gw.detector.InterferometerList(['H1', 'L1'])
interferometers.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] - 6,
)
for ifo in interferometers:
    ifo.minimum_frequency = minimum_frequency
    ifo.maximum_frequency = maximum_frequency
# SEBONRe waveform to inject
t, seobnre_waveform_time_domain = wf.seobnre_bbh_with_spin_and_eccentricity(
    parameters=injection_parameters,
    sampling_frequency=sampling_frequency,
    minimum_frequency=minimum_frequency - 10,
    maximum_frequency=maximum_frequency,
)
seobnre_wf_td, seobnre_wf_fd, max_overlap, index_shift, phase_shift = ovlp.maximise_overlap(
    seobnre_waveform_time_domain,
    comparison_waveform_frequency_domain,
    sampling_frequency,
    interferometers[0].frequency_array,
    interferometers[0].power_spectral_density,
)
seobnre_wf_fd = ovlp.zero_pad_frequency_domain_signal(
    seobnre_wf_fd, interferometers
)
# Inject the signal
interferometers.inject_signal(
    parameters=injection_parameters, injection_polarizations=seobnre_wf_fd
)
# Set up priors
priors = bb.gw.prior.BBHPriorDict()
for key in ['mass_1', 'mass_2', 'a_1', 'a_2']:
    priors.pop(key)
for key in ['theta_jn', 'theta_jl', 'tilt_1', 'tilt_2']:
    priors[key] = 0
priors['mass_ratio'] = bb.core.prior.Uniform(name='mass_ratio', minimum=0.125, maximum=1, boundary='reflective')
priors['chirp_mass'] = bb.core.prior.Uniform(name='chirp_mass', minimum=9.0, maximum=69.9, unit='$M_{\\odot}$', boundary='reflective')
# eccentricity = LogUniform(name='eccentricity', minimum=1e-4, maximum=0.2, boundary='reflective')
priors['chi_1'] = bb.gw.prior.AlignedSpin(a_prior=bb.gw.prior.Uniform(0, 0.6), z_prior=bb.gw.prior.Uniform(-1, 1), name='chi_1', latex_label='$\\chi_1$', boundary='reflective')
priors['chi_2'] = bb.gw.prior.AlignedSpin(a_prior=bb.gw.prior.Uniform(0, 0.6), z_prior=bb.gw.prior.Uniform(-1, 1), name='chi_2', latex_label='$\\chi_2$', boundary='reflective')
priors['luminosity_distance'] = bb.core.prior.PowerLaw(alpha=2, name='luminosity_distance', minimum=50, maximum=1000, unit='Mpc', latex_label='$d_L$')
priors['dec'] = bb.core.prior.Cosine(name='dec', boundary='reflective')
priors['ra'] = bb.core.prior.Uniform(name='ra', minimum=0, maximum=2 * np.pi, boundary='periodic')
priors['theta_jn'] = bb.core.prior.Sine(name='theta_jn', boundary='reflective')
priors['psi'] = bb.core.prior.Uniform(name='psi', minimum=0, maximum=np.pi, boundary='periodic')
priors['phase'] = bb.core.prior.Uniform(name='phase', minimum=0, maximum=2 * np.pi, boundary='periodic')
priors['geocent_time'] = bb.core.prior.Uniform(
    name='phase', minimum=injection_parameters['geocent_time'] - (deltaT / 2),
    maximum=injection_parameters['geocent_time'] + (deltaT / 2)
)
# Likelihood
likelihood = bb.gw.likelihood.GravitationalWaveTransient(
    interferometers, comparison_waveform_generator, time_marginalization=True,
    phase_marginalization=True, distance_marginalization=True, priors=priors
)
# Launch sampler
result = bb.core.sampler.run_sampler(
    likelihood,
    priors,
    sampler="dynesty",
    npoints=2000,
    walks=200,
    outdir=outdir,
    label=label,
)
# Plot corner
result.samples = bb.gw.conversion.generate_posterior_samples_from_marginalized_likelihood(result.samples, likelihood)
result.plot_corner()
