### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 7c4adbf2-309d-44ed-9206-6d4c1470b635
begin
	using Plots
	
	x = range(0, 20π, length=1000)  
	wrapped_phase = mod.(x .+ π, 2π) .- π  
	unwrapped_phase = x  
	
	plot(x, unwrapped_phase, label="Unwrapped Phase", lw=4, xlabel="Time (ms)", ylabel="Phase (radians)", legend=:top, fontsize=14, labelfont=14, legendfont=14, tickfont=14, xguidefontsize=14, yguidefontsize=14)
	
	plot!(x, wrapped_phase, label="Wrapped Phase", lw=4, color=:red, fontsize=14, labelfont=14, legendfont=14, tickfont=14, xguidefontsize=14, yguidefontsize=14, grid=false, title="Unwrapped vs Wrapped Phases")
end

# ╔═╡ 5ff814bb-706c-4266-9aa5-cc5fe7486a5f
using StatsBase

# ╔═╡ 86886228-9e3a-4854-ae95-d5a0bf0442b7
using CSV, DataFrames, FFTW, Random

# ╔═╡ 60ae3c91-2764-46f0-b916-7474b9e466f0
using SpecialFunctions

# ╔═╡ 3b71e1f0-0c4a-11f1-0aa6-43a57fc47120
md"""
# **Data-driven analysis of seismic phase using circular statistics**
"""

# ╔═╡ 67f35fa0-dfa6-424d-ba6b-7e9005936ff1
md"""
Akshika Rohatgi, Andrey Bakulin, and Sergey Fomel
"""


# ╔═╡ f40c721d-f78c-41e8-97ed-efd0625d88c3
md"""
# **Abstract**
"""

# ╔═╡ a6d2d852-80d0-4074-83e2-53854ad1e5de
md"""
Recognizing seismic phase as a primary attribute in seismic processing workflows, we apply circular statistics, a robust data-driven approach for correcting phase distortions in prestack seismic data. Unlike traditional linear methods that struggle with wrapped phase and often defer phase diagnostics to the final processing stages, the proposed approach treats phase as a circular variable.

We compute the circular mean, variance, and von Mises concentration parameter directly from phase ensembles in the frequency domain. These parameters provide insights into phase stability and coherence without needing phase unwrapping or wavelet assumptions.

Synthetic tests using additive and multiplicative noise models confirm that phase distributions follow the von Mises distribution—an analog of the normal distribution for circular variables—with circular statistics reliably tracking the true phase even in low-signal-quality scenarios.

Field examples demonstrate how this framework can map phase behavior across frequency and offset, enabling the detection of coherence bands and assessing the impact of each processing step on phase fidelity. The proposed approach can be particularly valuable in land acquisition, where prestack data often exhibit a low signal-to-noise ratio. Circular statistics allow us to evaluate phase integrity at each frequency, facilitating novel data conditioning and acquisition design strategies.
"""


# ╔═╡ a400b618-d29d-46bd-87e0-09f4b92bc574
md"""
# **Introduction**
"""

# ╔═╡ cd96e0e5-e393-4f43-a7a2-ec0e15d95784
md"""
The phase component of seismic signals contains critical information for accurate subsurface imaging and for identifying structural and stratigraphic features (Oppenheim and Lim, 1981; Lichman, 1999; Ulrych et al., 2007; van der Baan and Fomel, 2009; Greer et al., 2019). In many applications—particularly migration, coherence analysis, and inversion—consistent phase alignment can be as important as, or even more important than, amplitude fidelity (Xie et al., 2016; Holt and Lubrano, 2020; Bakulin et al., 2022a, 2022b). As acquisition and processing move toward denser spatial sampling and larger data volumes, statistical analysis of phase variability becomes both feasible and increasingly necessary.

Conventional processing methods, including noise attenuation, surface-consistent scaling, and deconvolution (Taner and Koehler, 1981; Cary and Lorentz, 1993; Chan and Stewart, 1994; Meunier, 1999; Liu et al., 2006), often do not fully account for the true statistical nature of seismic phase, particularly under low signal-to-noise ratio (S/N) conditions or significant wavelet perturbations. A fundamental challenge is that phase is a circular variable bounded within [−π, π], violating the assumptions underlying standard linear statistics. As a result, conventional measures such as linear means and variances can lead to misleading interpretations (Bakulin et al., 2024).

This paper introduces a robust analytical framework based on circular statistics (Mardia and Jupp, 1999) to characterize seismic phase distributions across frequency. Treating phase as a circular random variable allows us to quantify both its average behavior (circular mean) and its dispersion or concentration (via the von Mises concentration parameter, κ), which can be interpreted as a frequency-dependent analog of S/N. This shift from linear to circular statistical analysis provides a more meaningful description of phase variability, particularly for noisy land data.

Unlike traditional workflows that evaluate phase quality primarily in stacked or migrated images, the proposed framework enables diagnostics directly on raw prestack data. Phase integrity can therefore be monitored throughout each processing stage. Because the analysis is performed frequency by frequency, it offers a powerful tool for assessing how processing affects phase fidelity—information that amplitude-based or stack-based S/N metrics often fail to capture.

This capability is especially important in land acquisition, where prestack data commonly exhibit low S/N, making conventional quality measures unreliable (Bakulin et al., 2022a). By directly evaluating phase coherence across frequency, the framework complements frequency-dependent processing approaches (Retailleau et al., 2014; Bakulin et al., 2019) and supports emerging techniques such as time–frequency phase masking (Bakulin et al., 2023), which estimate signal phase through local stacking. Here, we extend these ideas by analyzing the full statistical distribution of phase to estimate signal phase and quantify noise-induced perturbations.

We propose that phase statistics should be treated as primary diagnostic attributes, monitored alongside conventional S/N throughout the processing workflow. Such diagnostics can inform acquisition design, bandwidth evaluation, and the development of phase-aware processing strategies. Tests on both synthetic and field data demonstrate that the proposed approach yields stable, interpretable estimates of phase behavior even under severe noise conditions.
"""


# ╔═╡ a911602a-79e2-4b77-9be7-cdf2462dda18
md"""
# **When phase misleads: A visual case for rethinking analysis**
"""

# ╔═╡ dc709536-5ebb-4d6e-b196-1040262c5f7a
md"""
We begin with a simple yet illustrative synthetic example that highlights the challenges of phase analysis in noisy environments. Figure 1 shows a seismic ensemble containing three flat reflectors, each contaminated by additive Gaussian noise. Two noise levels are considered: −2 dB and −10 dB S/N. The underlying signal is identical across all traces, so any phase variability arises solely from noise.

Figures 1c and 1d present histograms of the estimated phase at 40 Hz and 60 Hz on a linear scale. In the absence of noise, the phase would be constant across traces, collapsing the distribution into a narrow spike (dashed red line). Even moderate noise, however, spreads the phase distribution. As the noise level increases, the contamination becomes more severe.

At 40 Hz (Figure 1c), the distribution appears roughly symmetric, suggesting small random perturbations around the true phase. Closer inspection, however, reveals that the distribution is broader and less centered than expected. At 60 Hz (Figure 1d), the effect is more pronounced: the phase distribution becomes bimodal, splitting into two distinct peaks.

This behavior is not simply random scatter but reflects the limitations of linear statistical analysis applied to circular variables. Linear methods implicitly treat phase as an ordinary real-valued quantity, ignoring its periodic nature. As a result, they can produce misleading interpretations when phase wraps around the [−π, π] interval, especially under low S/N conditions.
"""


# ╔═╡ ff4898a3-41dd-4933-8c33-a34e569bb36b
md"""
# **From confusion to clarity: Circular statistics explained visually**
"""

# ╔═╡ 4e7f871c-4850-4960-afdb-6160d07c04e2
md"""
To interpret the distorted phase patterns observed earlier, we must reconsider how phase behaves. Phase is not an ordinary scalar quantity; it is periodic and circular. When phase exceeds π or falls below −π, it wraps around the unit circle. This property fundamentally changes how averages and variability should be computed. Figure 2 illustrates this effect: while the unwrapped phase evolves smoothly, the wrapped phase—confined to [−π, π]—exhibits abrupt jumps. These discontinuities are not physical features of the signal but artifacts of representation, yet they severely complicate the use of conventional linear statistics.

Phase unwrapping is widely used in applications such as radar, InSAR, and optical interferometry, where continuous phase fields are required to estimate displacement or topography (Ghiglia and Romero, 1994; Zebker and Lu, 1998). In seismic data, however, wavefields are often noisy and complex, causing unwrapping algorithms to become unstable or unreliable. Moreover, the wrapped phase itself is a valid bounded quantity and should not necessarily be forced into a linear framework. Instead, it naturally lends itself to circular analysis.

Circular statistics provide the appropriate mathematical framework. Rather than treating phase as a value on a line, each phase angle θ is represented as a direction on the unit circle. For a set of N phases ``(\theta_1, \theta_2, \ldots, \theta_N)``, we compute the cosine and sine averages:
```math
C = \frac{1}{N}\sum_{i=1}^{N} \cos \theta_i, \qquad
S = \frac{1}{N}\sum_{i=1}^{N} \sin \theta_i
```

The length of the resultant vector is
```math
R = \sqrt{C^2 + S^2}, \qquad R \ge 0
```

The circular mean direction is
```math
\bar{\theta} = \operatorname{atan2}(S,\, C)
```

This formulation avoids classic linear pitfalls—for example, the average of 1° and 359° is correctly computed as 0°, not 180°.

Beyond the mean, phase coherence is quantified using the mean resultant length
```math
\bar{R} = R, \qquad 0 \le \bar{R} \le 1
```

Values close to 1 indicate strong phase alignment, whereas values near 0 indicate a highly scattered distribution.

The circular variance is defined as
```math
V = 1 - \bar{R}
```

Perfect alignment yields ``V = 0``, while uniformly distributed phases yield ``V \approx 1``. Unlike linear variance, this measure remains meaningful for wrapped data.

A key advantage of circular statistics is that these quantities—the circular mean, resultant length, and circular variance—can be computed directly from wrapped phase values without phase unwrapping or signal modeling. They therefore provide simple, robust descriptors of phase behavior analogous to mean and standard deviation in linear statistics, but adapted to angular data.

Table 1 summarizes the principal attributes used in this study. While the circular mean, resultant length, and variance apply to any angular distribution, the concentration parameter κ is specific to the von Mises distribution. Because κ can be estimated directly from ``\bar{R}``, it serves as a useful derived metric when a von Mises model is appropriate. This relationship becomes central in the next section.
"""

# ╔═╡ d4af8fc2-e90a-4e63-9294-a40498d9a7c6
md"""
| **Attribute** | **Symbol** | **Range** | **Interpretation** |
|:---|:---:|:---:|:---|
| Circular mean | ``\bar{\theta}`` | ``[-\pi, \pi]`` | Central phase direction |
| Mean resultant length | ``\bar{R}`` | ``[0, 1]`` | Degree of phase alignment (clustering strength) |
| Circular variance | ``V = 1 - \bar{R}`` | ``[0, 1]`` | Phase dispersion or incoherence |
| Concentration (von Mises) | ``\kappa`` | ``[0, \infty)`` | Inverse dispersion (analogous to S/N in phase space) |
"""

# ╔═╡ e3214f03-f2a4-4961-8a19-8972e55d58ba
md"""
# **Modeling phase distributions with von Mises: The circular Gaussian**
"""

# ╔═╡ 05619bbf-9e1b-4955-a7e6-2e849a0932d4
md"""
The circular mean and variance (equations 2 and 4, respectively) give us powerful ways to summarize phase behavior, but what if we want a complete, interpretable model of how phase values are distributed? In linear statistics, the normal (Gaussian) distribution serves this role. For circular data such as phase, the equivalent is the von Mises distribution (Mardia and Jupp, 1999).

This distribution provides a natural model for how wrapped phase values behave in the presence of common types of seismic noise. Both additive Gaussian and multiplicative random noise, typical in land seismic data, lead to phase perturbations that closely follow von Mises distributions (Bakulin et al., 2022b, 2024; Rohatgi et al., 2024a).

The von Mises probability density function is:
```math
f(\theta;\, \bar{\theta},\, \kappa) = \frac{1}{2\pi I_0(\kappa)} \exp\!\left[\kappa \cos(\theta - \bar{\theta})\right]
```

where ``\theta`` is the wrapped phase angle, ``\bar{\theta}`` is the circular mean representing the central direction, ``\kappa`` is the concentration parameter analogous to the inverse of variance, and ``I_0(\kappa)`` is the modified Bessel function of the first kind of order 0.

This distribution wraps a bell-shaped curve around the unit circle. For large ``\kappa``, it peaks tightly around the mean. For small ``\kappa``, it flattens out into a near-uniform circle.

The von Mises distribution fits real seismic data well because it captures the kind of phase variability induced by realistic random noise, not through theoretical assumptions, but directly measurable in ensembles of traces. This allows us to go beyond just detecting noise — we can statistically model it. The distribution provides two key parameters that matter most for analysis:

- The circular mean ``\bar{\theta}`` represents the direction of dominant or average phase and estimates the underlying signal phase.
- The concentration parameter ``\kappa``, specific to the von Mises distribution, quantifies how tightly the phase vectors are clustered around the mean. It serves as a phase-based analog of S/N.

Importantly, ``\kappa`` can be estimated from ``\bar{R}`` using the following approximations (Fisher et al., 1993):
```math
\kappa = \begin{cases}
2\bar{R} + \bar{R}^3 + \dfrac{5\bar{R}^5}{6}, & \bar{R} < 0.53 \\[8pt]
-0.4 + 1.39\bar{R} + \dfrac{0.43}{1 - \bar{R}}, & 0.53 \le \bar{R} < 0.85 \\[8pt]
\dfrac{1}{3\bar{R} - 4\bar{R}^2 + \bar{R}^3}, & \bar{R} \ge 0.85
\end{cases}
```

"""

# ╔═╡ 24e2c578-835b-475d-8947-709d9fa65677
function wrap_phase(unwrapped_phase)
    return mod.(unwrapped_phase .+ π, 2π) .- π
end

# ╔═╡ 2b06eb0b-51f2-45e3-af2e-6558a70d9bac
begin
	using Distributions, CircStats

	kappa_values = 0.001:0.1:20
	circ_var_values = map(κ -> 1.0 - circ_r(wrap_phase.(rand(VonMises(0.0, κ), 100_001))), kappa_values)

	plot(kappa_values, circ_var_values, xlabel="κ", ylabel="Circular Variance", title="Circular Variance vs. κ", legend=false, linewidth=4, color=:black, grid=false, tickfontsize=12, xguidefontsize=14, yguidefontsize=14, titlefontsize=16, dpi=600, size=(700, 450))
end

# ╔═╡ efd83afd-bde3-4377-a1fb-f84a55307435
md"""
In short, low variance (``V`` close to 0) corresponds to high concentration (large ``\kappa``), indicating tightly clustered, coherent phases. High variance (``V`` near 1) implies low concentration (small ``\kappa``), as in a uniform phase distribution. This relationship is visualized in Figure 5.

Theoretically, if phase unwrapping were reliably possible, one could fit a traditional Gaussian to the unwrapped phases and estimate standard deviation ``\sigma`` in radians or degrees, giving familiar uncertainty metrics. Bakulin et al. (2024) derived numerical correspondences between ``\kappa`` and ``\sigma`` for such cases. However, unwrapping is rarely feasible or stable for noisy seismic data in practice.

Therefore, becoming fluent in interpreting wrapped-domain quantities like circular variance (``V``) and the von Mises concentration parameter (``\kappa``) is essential. Both describe phase dispersion, but differently. Circular variance is a general-purpose metric that can be computed for any phase distribution, whereas ``\kappa`` is specific to the von Mises model. ``V`` is often more robust and intuitive for practical diagnostics in field data due to its bounded range between 0 and 1. In contrast, ``\kappa`` ranges from 0 to ∞ and can exhibit large fluctuations between neighboring samples in noisy settings.

Nonetheless, when the von Mises model is assumed, ``V`` and ``\kappa`` are tightly linked through the mean resultant length ``\bar{R}`` and are analytically or graphically interconvertible (see Figure 5 and Table 1).

In short, the von Mises distribution gives us more than a fit — it offers a framework. It allows us to interpret phase coherence, model uncertainty, and detect reliable signal content across frequencies. With this model in hand, we are now ready to return to our earlier examples, this time through the lens of circular statistics, and see what was hidden in plain sight.
"""

# ╔═╡ 05fd2bdc-299a-4470-baff-e9f60887bb1d
function variance_to_kappa(V::Real)
    R = 1 - V
    if R < 1e-8
        return 0.0
    elseif R < 0.53
        return 2R + R^3 + (5R^5)/6
    elseif R < 0.85
        return -0.4 + 1.39R + 0.43/(1 - R)
    else
        return 1 / (R^3 - 4R^2 + 3R)
    end
end

# ╔═╡ 4977c317-74cc-43ef-9f4b-ccbafbbc0496
function add_circle_axes!(p; r=1.0)
    tt = range(0, 2π; length=361)
    plot!(p, r*cos.(tt), r*sin.(tt), lc=:black, lw=3, label=false)
    plot!(p, [-r, r], [0, 0], lc=:black, lw=3, label=false)
    plot!(p, [0, 0], [-r, r], lc=:black, lw=3, label=false)
    return p
end

# ╔═╡ 99e07135-36d9-4820-bb89-29a4194715f6
function plot_vectors(θ; step=5, col=:lightgreen, ttl="")
	xs, ys = cos.(θ), sin.(θ)
	R_bar = circ_r(θ)
	mean_dir, = circ_mean(θ)
	R_bar_vector = (R_bar*cos(mean_dir), R_bar*sin(mean_dir))
	idx = 1:step:length(xs)

	p = quiver(zeros(length(idx)), zeros(length(idx)); quiver=(xs[idx], ys[idx]), aspect_ratio=1, legend=false, grid=false, xlims=(-1.1,1.1), ylims=(-1.1,1.1), arrow=true, lw=1, linecolor=col, axis=false, title=ttl)

	add_circle_axes!(p)
	annotate!(p, [(1.05, 0.00, text("0", 14)), (-1.08, 0.00, text("-π,\n π", 14, :center)), (0.00, 1.08, text("π/2", 14, :center)), (0.00, -1.10, text("-π/2", 14, :center))])

	quiver!(p, [0.0], [0.0]; quiver=([R_bar_vector[1]], [R_bar_vector[2]]), arrow=true, lw=6, linecolor=:green, label=false)
	return p
end

# ╔═╡ 21204036-213c-4a95-afed-b9e861238584
function c_hist(θ; nbins=30, colR=:green, ttl="")

    # Histogram edges on [-π, π]
    edges = range(-π, π; length=nbins+1)
    h = fit(Histogram, θ, edges)

    counts  = h.weights
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    widths  = diff(edges)

    # Normalize radii 
    maxc = maximum(counts)
    r = maxc > 0 ? counts ./ maxc : zeros(length(counts))
    
    # Base plot
    p = plot(aspect_ratio=1, legend=false, grid=false,
             xlims=(-1.1, 1.1), ylims=(-1.1, 1.1),
             title=ttl, axis=false)

    for (c, w, rr) in zip(centers, widths, r)
        t1 = c - w/2
        t2 = c + w/2
        tseg = range(t1, t2; length=12)

        x = vcat(0.0, rr*cos.(tseg), 0.0)
        y = vcat(0.0, rr*sin.(tseg), 0.0)

        plot!(p, x, y, seriestype=:shape, fc=:white, lc=:black, lw=1, label=false)
    end

    add_circle_axes!(p)

    annotate!(p, [(1.05, 0.00, text("0", 14)), (-1.08, 0.00, text("-π,\n π", 14, :center)), (0.00, 1.08, text("π/2", 14, :center)), (0.00, -1.10, text("-π/2", 14, :center))])
    
    R_bar = circ_r(θ)
	mean_dir, = circ_mean(θ)
	R_bar_vector = (R_bar*cos(mean_dir), R_bar*sin(mean_dir))
    # Resultant vector overlay
    quiver!(p, [0.0], [0.0]; quiver=([R_bar_vector[1]], [R_bar_vector[2]]), arrow=true, lw=6, linecolor=:green, label=false)

    return p
end

# ╔═╡ 3431bac7-0053-4597-977c-a3373bb07803
begin
	V1, V2 = 0.8, 0.1
	κ1, κ2 = variance_to_kappa(V1), variance_to_kappa(V2)
	R1 = 1-V1
	R2 = 1-V2
	
	d1, d2 = VonMises(0.0, κ1), VonMises(0.0, κ2)
	phase1 = rand(d1, 1000)
	phase2 = rand(d2, 1000)
end

# ╔═╡ 211309dc-b246-4182-ab11-d34a52321b5d
begin
	va = plot_vectors(phase1; ttl="R̄=$(round(R1,digits=2)), V=$V1")
	vb = plot_vectors(phase2; ttl="R̄=$(round(R2,digits=2)), V=$V2")
	ca = c_hist(phase1; nbins=30, ttl="R̄=$(round(R1,digits=2)), V=$V1")
	cb = c_hist(phase2; nbins=30, ttl="R̄=$(round(R2,digits=2)), V=$V2")

	plot(va, vb, ca, cb, layout=(2,2), size=(900,900))
end

# ╔═╡ e98e9a2d-aefe-4033-a6b6-17cfad1bdf6e
md"""
# **Application of circular statistics to the additive white noise example**
"""

# ╔═╡ 0169b5c9-0ae6-414a-9442-51b79a05b6ff
md"""
Figures 1c–1h and 1k–1p apply circular statistics to the synthetic additive white noise example, accumulating phase ensembles from 1000 traces and analyzing them frequency by frequency. The circular mean aligns precisely with the true signal phase, and the observed distributions fit well with the von Mises model (equation 5).

Traditional linear histograms (Figures 1c–1d) fail to capture true phase behavior. The apparent bimodality at 60 Hz (Figure 1d) is an artifact of linear representation, not a physical phenomenon — seismic phase is fundamentally circular and bounded within ``[-\pi, \pi]``.

Rose diagrams (Figures 1g–1h, 1o–1p) resolve this ambiguity: what appeared bimodal on a line forms a single coherent cluster around the true phase direction. Each unit vector represents a trace's phase at a given frequency; their average direction robustly estimates the dominant phase with no unwrapping or modeling required.

As expected from white noise, variance (``\kappa`` or ``V``) remains constant across frequencies — only the signal phase itself rotates naturally. Circular statistics thus reveal structure and coherence that traditional methods obscure.
"""

# ╔═╡ 419bc868-69ec-41b6-8d1e-8aaf3d5a29a5
T = 0.002

# ╔═╡ a0787960-a985-4a22-bfa2-fbeb86dda5c1
fs = 1/T

# ╔═╡ a6c85b6d-ad77-48da-896c-60f3db9ede1e
num_wiggles=100

# ╔═╡ d8337188-f09a-4155-8157-4154bed53f70
#clean traces
begin
    klauder = DataFrame(CSV.File("Klauder_wavelet.csv")).signal[300:800]
	clean_traces = repeat(klauder, 1, num_wiggles)
end

# ╔═╡ ad9e5ebf-991e-4141-a228-b938b1f37979
t = range(0, step=T, length=length(klauder))

# ╔═╡ 940cef5f-027c-43b4-8aec-5e74be24dd16
function add_noise(signal, num_wiggles, desired_SNR_dB)
    noise_power = var(signal) / 10^(desired_SNR_dB / 10)
    noisy_traces = zeros(Float64, length(signal), num_wiggles)
    for i in 1:num_wiggles
        noisy_traces[:,i] = signal .+ sqrt(noise_power) .* randn(size(signal))
    end
    return noisy_traces
end

# ╔═╡ 9ad3274c-c268-4a51-b004-12ae2e46bc6c
begin
	add_noise_traces2dB = add_noise(klauder, num_wiggles, -2)
	add_noise_traces10dB = add_noise(klauder, num_wiggles, -10)
end

# ╔═╡ 33a629b6-e4db-4dad-acdf-688eac1a904f
function trace_spectrum(data, dt; extras=false)
    nt, ntr = size(data)
    fft_trace = fft(data[:, 1])
    f_pos = fft_trace[2:end]
    freq = fftfreq(length(f_pos), 1/dt)
    pos_inds = findall(x -> x ≥ 0, freq)
    nfreq = length(pos_inds)
    
    phase_vals = zeros(Float32, nfreq, ntr)
    amplitude_vals = zeros(Float32, nfreq, ntr)
    f_zero = zeros(ComplexF64, 1, ntr)
    positive_freq = freq[pos_inds]

    for k in 1:ntr
        fft_trace = fft(data[:, k])
        f_zero[1, k] = fft_trace[1]
        f_pos = fft_trace[2:end]
        phase_vals[:, k] = angle.(f_pos[pos_inds])
        amplitude_vals[:, k] = abs.(f_pos[pos_inds])
    end

    return extras ? (phase_vals, amplitude_vals, positive_freq, f_zero) :    (phase_vals, amplitude_vals)
end

# ╔═╡ 8e7dc06e-2618-41ae-bfd8-26183a68a06b
begin
    phase_clean, amplitude_clean, positive_freq, f_zero = trace_spectrum(clean_traces, T, extras=true)
    phase_add_2dB,  amplitude_add_2dB  = trace_spectrum(add_noise_traces2dB,  T)
    phase_add_10dB, amplitude_add_10dB = trace_spectrum(add_noise_traces10dB, T)
end

# ╔═╡ bf54d506-2303-4a05-966f-6b18f0bb27bf
function phase_stats(phase_matrix)
    nfreq = size(phase_matrix, 1)
    mean_vals = zeros(Float64, nfreq)
    var_vals  = zeros(Float64, nfreq)
    for i in 1:nfreq
        row = phase_matrix[i, :]
        mean_vals[i], = circ_mean(row)
        var_vals[i]   = 1 - circ_r(row)
    end
    return mean_vals, var_vals
end

# ╔═╡ b6eb8779-a9fc-4bd0-9af1-df96f52db306
begin
	predicted_phase_2dB  = repeat(phase_stats(phase_add_2dB)[1],  1, num_wiggles)
	predicted_phase_10dB = repeat(phase_stats(phase_add_10dB)[1], 1, num_wiggles)
end

# ╔═╡ d613df3f-1ec1-465d-aa0f-ce71d2d05276
function von_mises_pdf(θ, μ, κ)
    1 / (2π * besseli(0, κ)) * exp(κ * cos(θ - μ))
end

# ╔═╡ 1f7fb95a-253a-4d9a-ba11-0d5795826cac
#phase distributions
function phase_dist(phases, clean_phases, predicted_phases, freq_idx, freq_label)
    mean_c, = circ_mean(phases[freq_idx, :])
	kappa_c = circ_kappa(phases[freq_idx, :])
    vm_dist = VonMises(mean_c, kappa_c)
    θ = range(-π, π, length=1000)
    pdf_vals = [von_mises_pdf(x, mean_c, kappa_c) for x in θ]
    ymax = maximum(pdf_vals)

    h = histogram(phases[freq_idx, :], fill="#6161f4", bins=100, label="$freq_label Hz", normalize=:true, xguidefontsize=12, yguidefontsize=12, grid=false)
    plot!(h, θ, pdf_vals, label="vonMises Fit", xlabel="θ (radians)", ylabel="Density", lw=4, color=:red)
    plot!(h, [predicted_phases[freq_idx,1], predicted_phases[freq_idx,1]], [0, ymax],
          color="#8ed973", lw=6, label="Predicted Phase", grid=false)
    plot!(h, [clean_phases[freq_idx,1], clean_phases[freq_idx,1]], [0, ymax],
          color=:red, lw=2, linestyle=:dash, label="True Phase", grid=false)
    return h
end


# ╔═╡ 52fa59b4-b583-477a-8a3c-d2555f464cba
function circ_phase_dis(phases, μ, κ, clean_phase, predicted_phase;
                        title_str="", fillcol="#6161f4", lw=3.5, nbins=50)
    edges = range(-π, π; length=nbins+1)
    h = fit(Histogram, phases, edges)
    counts = h.weights
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    widths = diff(edges)
    maxc = maximum(counts)
    r = maxc > 0 ? counts ./ maxc : zeros(length(counts))
    p = plot(aspect_ratio=1, legend=:topright, grid=false,
             xlims=(-1.1,1.1), ylims=(-1.1,1.1),
             title=title_str, axis=false)
    # histogram wedges
    for (c, w, rr) in zip(centers, widths, r)
        tseg = range(c - w/2, c + w/2; length=12)
        x = vcat(0.0, rr*cos.(tseg), 0.0)
        y = vcat(0.0, rr*sin.(tseg), 0.0)
        plot!(p, x, y,
              seriestype=:shape,
              fc=fillcol,
              lc=:black,
              lw=1,
              label=false)
    end
    add_circle_axes!(p)
    # von Mises fit
    θ = range(-π, π, length=2000)
    rv = @. exp(κ * cos(θ - μ)) / (2π * besseli(0, κ))
    rv = rv ./ maximum(rv)
    plot!(p, rv .* cos.(θ), rv .* sin.(θ),
          color=:red,
          lw=lw,
          label="VM fit")
    # predicted and true phase arrows
    for (angle, col, lab, ls) in [
        (predicted_phase, "#8ed973", "Predicted", :solid),
        (clean_phase, :red, "True phase", :dash)
    ]
        aw = mod(angle + π, 2π) - π
        plot!(p, [0, cos(aw)], [0, sin(aw)],
              lw=3,
              color=col,
              linestyle=ls,
              arrow=:arrow,
              label=lab)
    end
    annotate!(p, [
        ( 1.05,  0.00, text("0", 14)),
        (-1.08,  0.00, text("-π,\n π", 14, :center)),
        ( 0.00,  1.08, text("π/2", 14, :center)),
        ( 0.00, -1.10, text("-π/2", 14, :center))
    ])
    return p
end

# ╔═╡ 41eda86c-1f12-4439-ba96-c197ded4e8ba
function reconstruct_signal(amplitude, predicted_phases, f_dc, num_wiggles, signal_length)
    new_f = amplitude .* exp.(1.0im .* predicted_phases)
    
    mirror = hcat([reverse(conj(new_f[:, i])) for i in 1:num_wiggles]...)
    dc_row = transpose(fill(f_dc, num_wiggles))
    full_f  = vcat(dc_row, new_f, mirror)
    
    return hcat([real(ifft(full_f[:, i])) for i in 1:num_wiggles]...)
end

# ╔═╡ 30afee40-8618-4a9a-b4ca-cc186496c47b
begin
	reconstructed_2dB  = reconstruct_signal(amplitude_clean, predicted_phase_2dB,  f_zero[1], num_wiggles, size(add_noise_traces2dB, 1))
	reconstructed_10dB = reconstruct_signal(amplitude_clean, predicted_phase_10dB, f_zero[1], num_wiggles, size(add_noise_traces10dB, 1))
end

# ╔═╡ 66745e37-96a2-4fc6-bd8f-106c4b7d0c17
begin
	common = (legend=false, ylabel="Time (ms)", xlabel="Trace Index",
	          xguidefontsize=12, yguidefontsize=12, yflip=true)


	# Wiggle plots
	p1 = plot(; common..., title="Noisy Traces -2 dB")
	for i in 1:size(add_noise_traces2dB, 2)
	    plot!(p1, add_noise_traces2dB[:, i] .+ i, t, linecolor=:black)
	end

	p2 = plot(; common..., title="Corrected Traces")
	for i in 1:size(reconstructed_2dB, 2)
	    plot!(p2, reconstructed_2dB[:, i] .+ i, t, linecolor=:black)
	end

	# Phase heatmaps
	p3 = heatmap(1:num_wiggles, positive_freq, phase_add_2dB,
	             yflip=true, ylabel="Frequency (Hz)", xlabel="Trace",
	             cmap=:coolwarm, title="Noisy Phases", ylims=(0,100), cbar=false)

	p4 = heatmap(1:num_wiggles, positive_freq, predicted_phase_2dB,
	             yflip=true, ylabel="Frequency (Hz)", xlabel="Trace",
	             cmap=:coolwarm, title="Predicted Phases", ylims=(0,100), cbar=false)

	plot(p1, p2, 
	     p3, p4, 
	     layout=(2,2), size=(1000,700))
end

# ╔═╡ b3e3ff9f-a2c3-47c9-bb6b-6ba410680b59
begin
    h1 = phase_dist(phase_add_2dB, phase_clean, predicted_phase_2dB, 41, "40 Hz")
    h2 = phase_dist(phase_add_2dB, phase_clean, predicted_phase_2dB, 61, "60 Hz")

    c1 = circ_phase_dis(phase_add_2dB[41,:], circ_mean(phase_add_2dB[41,:])[1], circ_kappa(circ_r(phase_add_2dB[41,:])), phase_clean[41,1], predicted_phase_2dB[41,1]; title_str="40 Hz")

    c2 = circ_phase_dis(phase_add_2dB[61,:], circ_mean(phase_add_2dB[61,:])[1], circ_kappa(circ_r(phase_add_2dB[61,:])), phase_clean[61,1], predicted_phase_2dB[61,1]; title_str="60 Hz")

    plot(h1, h2, c1, c2; layout=(2,2), size=(1000,700), plot_title="Additive noise −2 dB")
end

# ╔═╡ 4759331e-8568-479a-810a-d3953ce56554
begin

	# Wiggle plots
	p5 = plot(; common..., title="Noisy Traces -10 dB")
	for i in 1:size(add_noise_traces10dB, 2)
	    plot!(p5, add_noise_traces10dB[:, i] .+ i, t, linecolor=:black)
	end

	p6 = plot(; common..., title="Phase-only Corrected Traces")
	for i in 1:size(reconstructed_10dB, 2)
	    plot!(p6, reconstructed_10dB[:, i] .+ i, t, linecolor=:black)
	end

	# Phase heatmaps
	p7 = heatmap(1:num_wiggles, positive_freq, phase_add_10dB,
	             yflip=true, ylabel="Frequency (Hz)", xlabel="Trace",
	             cmap=:coolwarm, title="Noisy Phases", ylims=(0,100), cbar=false)

	p8 = heatmap(1:num_wiggles, positive_freq, predicted_phase_10dB,
	             yflip=true, ylabel="Frequency (Hz)", xlabel="Trace",
	             cmap=:coolwarm, title="Predicted Phases", ylims=(0,100), cbar=false)

	plot(p5, p6, 
	     p7, p8, 
	     layout=(2,2), size=(1000,700))
end

# ╔═╡ 270901e6-eecd-47c0-94ef-05a440f1387b
begin
    h3 = phase_dist(phase_add_10dB, phase_clean, predicted_phase_10dB, 41, "40 Hz")
    h4 = phase_dist(phase_add_10dB, phase_clean, predicted_phase_10dB, 61, "60 Hz")

    c3 = circ_phase_dis(phase_add_10dB[41,:], circ_mean(phase_add_10dB[41,:])[1], circ_kappa(circ_r(phase_add_10dB[41,:])), phase_clean[41,1], predicted_phase_10dB[41,1]; title_str="40 Hz")

    c4 = circ_phase_dis(phase_add_10dB[61,:], circ_mean(phase_add_10dB[61,:])[1], circ_kappa(circ_r(phase_add_10dB[61,:])), phase_clean[61,1], predicted_phase_10dB[61,1]; title_str="60 Hz")

    plot(h2, h3, c3, c4; layout=(2,2), size=(1000,700), plot_title="Additive noise −2 dB")
end

# ╔═╡ 767aedf5-2c0e-4a77-955f-638e7d4f6764
md"""
# **Phase analysis under multiplicative noise**
"""

# ╔═╡ 96d56984-cca3-4d1e-aa21-5eb7ef3ab670
md"""
To test the robustness of circular statistics under more realistic conditions, we apply a multiplicative noise model introduced by Bakulin et al. (2022b, 2023), specifically designed to induce phase perturbations representative of scattering and near-surface het-erogeneity in land seismic data. Unlike linear models that define noise using standard deviation (σ), this approach uses a frequency-dependent  concentration  parameter  (κ)  to  directly  control  the  spread of phase distributions and better capture the non-uniform, frequency-dependent nature of seismic noise.

We invert κ from seismic trace ensembles using the mean resultant length –R, applying standard approximations from Fisher et  al.  (1993)  and  Mardia  and  Jupp  (1999)  (Figure  6a).  The  resulting κ( f ) curve quantifies phase alignment as a function of frequency, essentially mapping a phase-based analog to S/N (Figure 6b).

Conceptually, κ serves a role analogous to S/N. High κ values indicate stable, well-aligned phases and low phase noise. Low κ reflects scattered phases and high phase noise. This trend is evident  in  the  phase  histograms  at  20,  40,  and  60  Hz  (Figures 6d–6F); however, the behavior becomes far more trans-parent  and  unambiguous  in  the  corresponding  rose  plots  (Figures 6g–6i). As frequency increases, phase vectors become more dispersed and κ decreases, consistent with the design of the multiplicative noise model.
"""

# ╔═╡ b925dd01-27d3-4364-a9c6-6e15eab1731a
function multiplicative_noise(phase_clean, amp_clean, f_zero, kappa_in, num_wiggles; κmax=500.0)
    perturbed_phases = similar(phase_clean)

    for i in 1:size(phase_clean, 1)
        κ = kappa_in[i]
        κ = isfinite(κ) ? clamp(κ, 0.0, κmax) : κmax  # <-- key line

        vm = VonMises(0.0, κ)
        noise = rand(vm, size(phase_clean, 2))

        perturbed_phases[i, :] .= phase_clean[i, :] .+ noise
    end

    return reconstruct_signal(amp_clean, perturbed_phases, f_zero[1], num_wiggles, size(phase_clean, 1))
end

# ╔═╡ 7051b57d-cfc3-4b12-81c5-fae602f9090b
begin
    sig = 2π .* positive_freq * 0.004   ##### CHANGE ME #####
    kappa_in = zeros(length(sig))
    for (i, M) in enumerate(sig)
        phase = angle.(exp.(1.0im .* randn(length(sig)) .* M))
        kappa = circ_kappa(phase)
        kappa_in[i] = kappa
    end
	mul_noisy_traces = multiplicative_noise(phase_clean, amplitude_clean, f_zero, kappa_in, num_wiggles; κmax=500.0)
end

# ╔═╡ c585f48f-4522-4201-9cf8-5ccad1ab1aab
phase_mul,  amplitude_mul  = trace_spectrum(mul_noisy_traces,  T)

# ╔═╡ 3c0e4d1d-6cb5-4e41-a8b3-6ac72b52b5dd
predicted_phase_mul  = repeat(phase_stats(phase_mul)[1],  1, num_wiggles)

# ╔═╡ 109f17f7-4186-4e30-a635-aaac90263ccb
reconstructed_mul  = reconstruct_signal(amplitude_clean, predicted_phase_mul,  f_zero[1], num_wiggles, size(mul_noisy_traces, 1))

# ╔═╡ f1b88a07-20df-4aab-a359-1f577c3e5d74
begin
    # Multiplicative noisy traces
    p = plot(legend=false, xlabel="Trace Index", ylabel="Time (ms)", yflip=true)
    for i in 1:size(mul_noisy_traces, 2)
        plot!(p, mul_noisy_traces[:, i] .+ i, t, linecolor=:black, 
			  title="Multiplicative noise")
    end

    # Reconstructed traces
    pm = plot(legend=false, xlabel="Trace Index", ylabel="Time (ms)", yflip=true)
    for i in 1:size(reconstructed_mul, 2)
        plot!(pm, reconstructed_mul[:, i] .+ i, t, linecolor=:black,
			 title="Phase-only corrected traces")
    end

    # Subplot
    plot(p, pm, layout=(1,2), size=(900,400))
end

# ╔═╡ e0cf9298-8430-428c-86d5-237885138a46
md"""
At  20  Hz  (Figures  6d  and  6g),  phase  vectors  are  tightly  clustered, and the circular mean closely aligns with the true phase. At  60  Hz  (Figures  6f  and  6i),  phase  scatter  increases,  and  the  circular mean drifts, signaling reduced reliability in phase estima-tion under higher perturbations.This analysis also underscores the importance of ensemble size.  In  high-noise  conditions,  small  ensembles  may  yield  unstable phase estimates. Even if linear plots show significant deviations  from  the  true  phase,  the  circular  representation  (e.g., Figure 6i) often reveals that the mean direction is statisti-cally meaningful. However, achieving convergence requires a sufficiently large number of traces. In an earlier paper (Rohatgi et al., 2024b), we proposed estimating the ensemble size using the  standard  deviation  of  the  residual  phase,  which  requires  knowing the true phase. In contrast, the approach proposed here  infers  phase  stability  and  confidence  directly  from  raw  phase  distributions,  without  the  need  for  phase  unwrapping  or prior wavelet information.

In summary, circular statistics quantify phase coherence and provide  practical  guidance  for  designing  ensembles  that  yield  robust phase estimates, even in extremely noisy seismic data.
"""

# ╔═╡ d187926c-47f8-46e2-9767-e88bc4b5f121
begin
    h5 = phase_dist(phase_mul, phase_clean, predicted_phase_mul, 21, "20 Hz")
    h6 = phase_dist(phase_mul, phase_clean, predicted_phase_mul, 41, "40 Hz")
    h7 = phase_dist(phase_mul, phase_clean, predicted_phase_mul, 61, "60 Hz")


    plot(h5, h6, h7; layout=(1,3),size=(1300,400), plot_title="Multiplicative noise", ylims=(0,1))
end

# ╔═╡ 94e90594-aa3f-41be-becc-bfdbec42ff6c
begin
        c5 = circ_phase_dis(phase_mul[21,:], circ_mean(phase_mul[21,:])[1], circ_kappa(circ_r(phase_mul[21,:])), phase_clean[21,1], predicted_phase_mul[21,1]; title_str="20 Hz")

        c6 = circ_phase_dis(phase_mul[41,:], circ_mean(phase_mul[41,:])[1], circ_kappa(circ_r(phase_mul[41,:])), phase_clean[41,1], predicted_phase_mul[41,1]; title_str="40 Hz")

        c7 = circ_phase_dis(phase_mul[61,:], circ_mean(phase_mul[61,:])[1], circ_kappa(circ_r(phase_mul[61,:])), phase_clean[61,1], predicted_phase_mul[61,1]; title_str="60 Hz")

    plot(c5, c6, c7; layout=(1,3), legend=false, size=(1200,400))
end

# ╔═╡ e5bedb7e-3cb2-48a5-b27c-196102bfcfbb
md"""
# **Conclusions**
"""

# ╔═╡ d190ddf1-b320-4ba9-8da1-254e682259f7
md"""
We have presented circular statistics as a practical and statisti-cally rigorous framework for seismic phase analysis. By treating phase as a directional variable confined to the interval [–π, π], we avoid  the  distortions  and  ambiguities  inherent  in  conventional  linear analysis, particularly under low signal-to-noise conditions.

The von Mises distribution provides a natural model for phase variability in noisy seismic data, enabling direct quantification of central tendency (via the circular mean) and coherence (via circular variance or κ) from raw phase values without the need for unwrap-ping or wavelet assumptions. The circular mean can reliably track the underlying signal phase for practically important additive and multiplicative  random  noise  cases,  even  in  highly  scattered  or  noise-dominated  environments.  While  similar  phase  estimates  could be obtained through local stacking, the circular statistics approach offers a deeper advantage: it recovers and visualizes the full absolute phase distribution, frequency by frequency, delivering both the mean phase and the associated phase spread. This richer quantitative insight into phase stability allows direct assessment of phase quality, which is critical for designing data-driven work-flows,  optimizing  ensemble  size,  and  enhancing  frequency-dependent processing strategies.

In comparison to traditional linear phase diagnostics, circular statistics,  as  visualized  through  rose  plots,  provide  clearer  and  more interpretable insights into phase behavior across frequency, offset, and processing parameters. All key quantities — circular mean, variance, and concentration — can be computed directly from the phase ensembles, eliminating the need for unwrapping, wavelet assumptions, or manual tuning. This reduces subjective bias and enables reproducible, objective data-driven analysis.

Phase  coherence  varies  with  frequency  across  the  seismic  spectrum,  presenting  opportunities  for  frequency-dependent  processing.  These  strategies  are  increasingly  critical  in  modern  broadband  workflows  to  enhance  resolution  and  recover  subtle  structural details. Whether used for quality control, preprocessing, or  as  a  foundation  for  advanced  imaging  techniques,  circular  statistics can serve as a robust and interpretable basis for phase-aware seismic analysis.
"""

# ╔═╡ 423e4c92-8645-4174-bec6-a568921b0d72
md"""
# **Acknowledgments**
"""

# ╔═╡ 89bbd044-b74b-4245-ba0b-20b5673d15f5
md"""
We  acknowledge  Fairfield  Geotechnologies  for  granting  permission to use the data presented in this study.
"""

# ╔═╡ 237f7c93-a97f-4e35-8071-1964bc78c207
md"""
# **References**
"""

# ╔═╡ 745fb91a-00f7-4f89-976f-b508dce0238b
md"""
Bakulin, A., D. Neklyudov, and I. Silvestrov, 2024, The impact of receiver arrays on suppressing seismic speckle scattering noise caused by the meter-scale near-surface heterogeneity: Geophysics, 89, no. 6, V551–V561, https://doi.org/10.1190/geo2023-0489.1.

Bakulin, A., D. Neklyudov, and I. Silvestrov, 2023, Seismic time-fre-quency masking for suppression of seismic speckle noise: Geophysics, 
88, no. 5, V371–V385, https://doi.org/10.1190/geo2022-0779.1.Bakulin, A., I. Silvestrov, and M. Protasov, 2022a, Signal-to-noise ratio computation  for  challenging  land  single-sensor  seismic  data:  Geophysical  Prospecting,  70,  no.  3,  629– 638,  https://doi.org/10.1111/1365-2478.13178.

Bakulin,  A.,  D.  Neklyudov,  and  I.  Silvestrov,  2022b,  Multiplicative  random seismic noise caused by small-scale near-surface scattering and its transformation during stacking: Geophysics, 87, no. 5, V419–V435, https://doi.org/10.1190/geo2021-0830.1.

Bakulin, A., I. Silvestrov, M. Dmitriev, D. Neklyudov, M. Protasov, K. Gadylshin,  and  A.  Dolgov,  2020a,  Nonlinear  beamforming  for  enhancement of 3D prestack land seismic data: Geophysics, 85, no. 3, V283-V296, https://doi.org/10.1190/geo2019-0341.1.

Bakulin,  A.,  D.  Neklyudov,  and  I.  Silvestrov,  2020b,  Prestack  data  enhancement with phase corrections in time‐frequency domain guided by  local  multidimensional  stacking:  Geophysical  Prospecting,  68, no. 6, 1811–1818, https://doi.org/10.1111/1365-2478.12956.

Bakulin, A., I. Silvestrov, and M. Dmitriev, 2019, Adaptive multiscale processing  of  challenging  3D  seismic  data  for  first-break  picking,  FWI  and  imaging:  89th  Annual  International  Meeting,  SEG,  Expanded Abstracts, 3979–3984, https://doi.org/10.1190/segam2019-3214616.1.

Cary, P. W., and G. A. Lorentz, 1993, Four-component surface-consistent deconvolution:  Geophysics,  58,  no.  3,  383–392,  https://doi.org/10.1190/1.1443421.

Chan,  W.  K.,  and  R.  R.  Stewart,  1994,  3-D  f-k  filtering:  CREWES  Research Report, 6, 15.1–15.7.

Fisher, N. I., T. Lewis, and B. J. J. Embleton, 1993, Statistical analysis of spherical data: Cambridge University Press.

Ghiglia,  D.  C.,  and  L.  A.  Romero,  1994,  Robust  two-dimensional  weighted and unweighted phase unwrapping that uses fast transforms and iterative methods: Journal of the Optical Society of America. A, Optics, Image Science, and Vision, 11, no. 1, 107–117, https://doi.org/10.1364/JOSAA.11.000107.

Greer, S., S. Fomel, and M. Fry, 2019, Prestack phase corrections using local  seismic  attributes:  89th  Annual  International  Meeting,  SEG,  3949–3953, https://doi.org/10.1190/segam2019-3214515.1.

Holt, R., and A. Lubrano, 2020, Stabilizing the phase of onshore 3D seismic  data:  Geophysics,  85,  no.  6,  V473–V479,  https://doi.org/10.1190/geo2019-0695.1.

Lichman,  E.,  1999,  Automated  phase-based  moveout  correction:  69th Annual  International  Meeting,  SEG,  Expanded  Abstracts,  1150–1153, https://doi.org/10.1190/1.1820706.

Liu, C., Y. Liu, B. Yang, D. Wang, and J. Sun, 2006, A 2D multistage median filter to reduce random seismic noise: Geophysics, 71, no. 5, V105–V110, https://doi.org/10.1190/1.2236003.

Mardia, K. V., and P. E. Jupp, 1999, Directional statistics: John Wiley & Sons, https://doi.org/10.1002/9780470316979.

Meunier, J., 1999, 3D geometry, velocity filtering and scattered noise: 69th  Annual  International  Meeting,  SEG,  Expanded  Abstracts,  1216–1219, https://doi.org/10.1190/1.1820725.

Oppenheim,  A.  V.,  and  J.  S.  Lim,  1981,  The  importance  of  phase  in  signals:  Proceedings  of  the  IEEE,  69,  no.  5,  529–541,  https://doi.org/10.1109/PROC.1981.12022.

Retailleau,  M.,  R.  El  Asrag,  and  J.  Shorter,  2014,  Processing  land  broadband data: Challenges that Oman surveys present and how they are addressed: EAGE/SPG Workshop on Broadband Seismic, https://doi.org/10.3997/2214-4609.20141699 l.

Rohatgi, A., A. Bakulin, and S. Fomel, 2024a, Analyzing the impact of additive  and  multiplicative  noise  on  seismic  data  analysis:  Fourth  International  Meeting  for  Applied  Geoscience  &  Energy,  SEG/APPG,  Expanded  Abstracts,  2132–2136,  https://doi.org/10.1190/image2024-4086176.1.

Rohatgi, A., A. Bakulin, and S. Fomel, 2024b, Phase pilot recovery: A foundation  for  mitigating  speckle  scattering  noise  in  seismic  data:  Fourth  International  Meeting  for  Applied  Geoscience  &  Energy,  SEG/A APG,  Expanded  Abstracts,  2314 –2318,  https://doi.org/10.1190/image2024-4091774.1.

Taner,  M.  T.,  and  F.  Koehler,  1981,  Surface  consistent  corrections:  Geophysics, 46, no. 1, 17–22, https://doi.org/10.1190/1.1441133.

Ulrych, T. J., S. T. Kaplan, M. D. Sacchi, and E. Galloway, 2007, The essence of phase in seismic data processing and inversion: 77th Annual International Meeting, SEG, Expanded Abstracts, 1765–1769, https://doi.org/10.1190/1.2792834.

van der Baan, M., and S. Fomel, 2009, Nonstationary phase estimation using regularized local kurtosis maximization: Geophysics, 74, no. 6, A75–A80, https://doi.org/10.1190/1.3213533.

Xie, X.-B., H. Ning, and B. Chen, 2016, How scatterings from small-scale near-surface heterogeneities affecting seismic data and the quality of depth  image,  analysis  based  on  seismic  resolution:  76th  Annual  International  Meeting,  SEG,  Expanded  Abstracts,  4278–4282,  https://doi.org/10.1190/segam2016-13780223.1.

Zebker, H. A., and Y. Lu, 1998, Phase unwrapping algorithms for radar interferometry: Residue-cut, least-squares, and synthesis algorithms: Journal of the Optical Society of America. A, Optics, Image Science, and Vision, 15, no. 3, 586 –598, https://doi.org/10.1364/JOSAA.15.000586.

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CircStats = "2f6764a1-d620-4564-9394-76eb7c776766"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.10.15"
CircStats = "~1.0.4"
DataFrames = "~1.8.1"
Distributions = "~0.25.120"
FFTW = "~1.10.0"
Plots = "~1.41.4"
SpecialFunctions = "~2.6.1"
StatsBase = "~0.34.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "5dcc9740aa8b76b4d60da86e043fdd32b94a5cce"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "856ecd7cebb68e5fc87abecd2326ad59f0f911f3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.43"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CircStats]]
deps = ["Distributions", "HypothesisTests", "LinearAlgebra", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "ecfe2e9a260c4723026b4a71460cf0420def9e40"
uuid = "2f6764a1-d620-4564-9394-76eb7c776766"
version = "1.0.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "8010b6bb3388abe68d95743dcbea77650bb2eddf"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.3"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3e6d038b77f22791b8e3472b7c633acea1ecac06"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.120"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "173e4d8f14230a7523ae11b9a3fa9edb3e0efd78"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.14.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ee0585b62671ce88e48d3409733230b401c9775c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.22"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "7dd7173f7129a1b6f84e0f03e0890cd1189b0659"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.22+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Printf", "Random", "Roots", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "15319d5a767eb386bc4b702d5e025a0be62be293"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.11.5"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c9cbeda6aceffc52d8a0017e71db27c7a7c0beaf"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.5+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "063ef757a1e0e15af77bbe92be92da672793fd4e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.4"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "c5a07210bd060d6a8491b0ccdee2fa0235fc00bf"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.1.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "34f7e5d2861083ec7596af8b8c092531facf2192"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "8f528b0851b5b7025032818eb5abbeb8a736f853"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "8a433b1ede5e9be9a7ba5b1cc6698daa8d718f1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "8e45cecc66f3b42633b8ce14d431e8e57a3e242e"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a3c1536470bf8c5e02096ad4853606d7c8f62721"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.2"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "9297459be9e338e546f5c4bedb59b3b5674da7f1"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "6ab498eaf50e0495f89e7a5b582816e2efb95f64"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.54+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─3b71e1f0-0c4a-11f1-0aa6-43a57fc47120
# ╟─67f35fa0-dfa6-424d-ba6b-7e9005936ff1
# ╟─f40c721d-f78c-41e8-97ed-efd0625d88c3
# ╟─a6d2d852-80d0-4074-83e2-53854ad1e5de
# ╟─a400b618-d29d-46bd-87e0-09f4b92bc574
# ╟─cd96e0e5-e393-4f43-a7a2-ec0e15d95784
# ╟─a911602a-79e2-4b77-9be7-cdf2462dda18
# ╟─dc709536-5ebb-4d6e-b196-1040262c5f7a
# ╠═7c4adbf2-309d-44ed-9206-6d4c1470b635
# ╟─ff4898a3-41dd-4933-8c33-a34e569bb36b
# ╟─4e7f871c-4850-4960-afdb-6160d07c04e2
# ╟─d4af8fc2-e90a-4e63-9294-a40498d9a7c6
# ╟─e3214f03-f2a4-4961-8a19-8972e55d58ba
# ╟─05619bbf-9e1b-4955-a7e6-2e849a0932d4
# ╠═24e2c578-835b-475d-8947-709d9fa65677
# ╠═2b06eb0b-51f2-45e3-af2e-6558a70d9bac
# ╟─efd83afd-bde3-4377-a1fb-f84a55307435
# ╠═05fd2bdc-299a-4470-baff-e9f60887bb1d
# ╠═4977c317-74cc-43ef-9f4b-ccbafbbc0496
# ╠═99e07135-36d9-4820-bb89-29a4194715f6
# ╠═21204036-213c-4a95-afed-b9e861238584
# ╠═3431bac7-0053-4597-977c-a3373bb07803
# ╠═5ff814bb-706c-4266-9aa5-cc5fe7486a5f
# ╠═211309dc-b246-4182-ab11-d34a52321b5d
# ╟─e98e9a2d-aefe-4033-a6b6-17cfad1bdf6e
# ╟─0169b5c9-0ae6-414a-9442-51b79a05b6ff
# ╠═86886228-9e3a-4854-ae95-d5a0bf0442b7
# ╠═419bc868-69ec-41b6-8d1e-8aaf3d5a29a5
# ╠═a0787960-a985-4a22-bfa2-fbeb86dda5c1
# ╠═a6c85b6d-ad77-48da-896c-60f3db9ede1e
# ╠═d8337188-f09a-4155-8157-4154bed53f70
# ╠═ad9e5ebf-991e-4141-a228-b938b1f37979
# ╠═940cef5f-027c-43b4-8aec-5e74be24dd16
# ╠═9ad3274c-c268-4a51-b004-12ae2e46bc6c
# ╠═33a629b6-e4db-4dad-acdf-688eac1a904f
# ╠═8e7dc06e-2618-41ae-bfd8-26183a68a06b
# ╠═bf54d506-2303-4a05-966f-6b18f0bb27bf
# ╠═b6eb8779-a9fc-4bd0-9af1-df96f52db306
# ╠═60ae3c91-2764-46f0-b916-7474b9e466f0
# ╠═d613df3f-1ec1-465d-aa0f-ce71d2d05276
# ╠═1f7fb95a-253a-4d9a-ba11-0d5795826cac
# ╠═52fa59b4-b583-477a-8a3c-d2555f464cba
# ╠═41eda86c-1f12-4439-ba96-c197ded4e8ba
# ╠═30afee40-8618-4a9a-b4ca-cc186496c47b
# ╠═66745e37-96a2-4fc6-bd8f-106c4b7d0c17
# ╠═b3e3ff9f-a2c3-47c9-bb6b-6ba410680b59
# ╠═4759331e-8568-479a-810a-d3953ce56554
# ╠═270901e6-eecd-47c0-94ef-05a440f1387b
# ╟─767aedf5-2c0e-4a77-955f-638e7d4f6764
# ╟─96d56984-cca3-4d1e-aa21-5eb7ef3ab670
# ╠═b925dd01-27d3-4364-a9c6-6e15eab1731a
# ╠═7051b57d-cfc3-4b12-81c5-fae602f9090b
# ╠═c585f48f-4522-4201-9cf8-5ccad1ab1aab
# ╠═3c0e4d1d-6cb5-4e41-a8b3-6ac72b52b5dd
# ╠═109f17f7-4186-4e30-a635-aaac90263ccb
# ╠═f1b88a07-20df-4aab-a359-1f577c3e5d74
# ╟─e0cf9298-8430-428c-86d5-237885138a46
# ╠═d187926c-47f8-46e2-9767-e88bc4b5f121
# ╠═94e90594-aa3f-41be-becc-bfdbec42ff6c
# ╟─e5bedb7e-3cb2-48a5-b27c-196102bfcfbb
# ╟─d190ddf1-b320-4ba9-8da1-254e682259f7
# ╟─423e4c92-8645-4174-bec6-a568921b0d72
# ╟─89bbd044-b74b-4245-ba0b-20b5673d15f5
# ╟─237f7c93-a97f-4e35-8071-1964bc78c207
# ╟─745fb91a-00f7-4f89-976f-b508dce0238b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
