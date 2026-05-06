### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ d6f47eee-dddf-4454-ac98-7fdf6f25ff99
using CSV, DataFrames, FFTW, Random, Plots

# ╔═╡ a955bc18-5858-438f-b3d5-5821eceeb55a
using CircStats, StatsBase

# ╔═╡ a0cb1efc-57b9-4fa2-aead-3633cb6405f3
using SpecialFunctions

# ╔═╡ df2e462c-20dc-43a7-9b81-1be8245fcc96
using Distributions

# ╔═╡ 27413a40-36df-11f1-3735-f9f362703981
md"""
# Amplitude-Invariant Phase Masking for Coherence Recovery in Scattered Wavefields
"""

# ╔═╡ fd269c66-9602-4f90-81bb-ebb512c8bb35
md"""
**Akshika Rohatgi¹, Andrey Bakulin¹, and Sergey Fomel¹**

¹ Bureau of Economic Geology, University of Texas at Austin
"""

# ╔═╡ 8f3256ed-c1b2-4e84-9763-af9ce49677d0
md"""
Coherent summation of multichannel recordings relies on phase stability across measurements. In scattered wavefields, near-surface heterogeneities introduce record-dependent phase distortions that decorrelate otherwise coherent signals and degrade summation quality. Conventional approaches to recovering phase coherence implicitly tie phase estimation to signal amplitude, introducing bias where amplitude variability is large. We present a framework based on circular statistics that separates phase estimation from amplitude entirely, ensuring all records contribute equally to the representative phase regardless of their energy. Synthetic experiments confirm that the proposed method outperforms conventional approaches when amplitude variability is present, offering a practical path toward more robust coherence recovery in seismic imaging and other wave-based sensing applications.
"""

# ╔═╡ e0caa217-e6f3-4013-bf62-b4f82e392c2c
md""" 
## 1. Introduction

Phase is a fundamental carrier of information in signals¹. When phase remains stable across multichannel measurements, coherent summation enhances signal energy through constructive interference². The need to isolate and denoise phase has been recognized across domains ranging from speech processing³˒⁴ to InSAR deformation mapping⁵. More recently, the same need has been identified in multichannel seismic data, where near-surface scattering introduces strong phase distortions across recorded signals⁶.
 
In multichannel seismic imaging, signals are broadband and coherence must be evaluated across frequency rather than under single-frequency assumptions. Phase distortions introduced by scattering vary with frequency and receiver position, complicating coherent summation. Beamforming and stacking methods exploit phase alignment across receiver arrays to reinforce coherent arrivals. Recent advances have shown that selectively suppressing incoherent time–frequency components through phase masking can further improve coherence recovery⁷. However, existing masking approaches derive coherence weights from locally stacked or beamformed traces, implicitly coupling phase estimation to amplitude information. When amplitudes vary significantly across receivers, this coupling biases the phase estimate and limits coherence recovery precisely where it is most needed.
 
In this letter, we introduce an amplitude-invariant phase masking framework grounded in circular statistics. The method decouples phase estimation from amplitude entirely. Building on the work of Bakulin et al.⁷, the proposed approach applies a unit-magnitude operator in the complex frequency domain: at each frequency, the phase of each record is replaced with a representative estimate derived from neighboring measurements using the circular mean. We demonstrate that this framework stabilizes phase in scattered wavefields and improves coherence in the presence of both phase perturbations and amplitude variability.
"""

# ╔═╡ 96f1b4dd-3e4c-493a-8b0e-57d11d7dea4a
begin
	T = 0.002
	fs = 1/T
	num_wiggles= 5000
end

# ╔═╡ 4153b364-dd01-4ddf-b556-4a2a592049d7
begin
using HTTP
file_id = "1suyqkO5jKfG3As5iub7iuWMDFXibwLH4"
url = "https://drive.google.com/uc?export=download&id=$file_id"
response = HTTP.get(url)
klauder = DataFrame(CSV.File(IOBuffer(response.body))).signal[500:672]
clean_traces = repeat(klauder, 1, num_wiggles);
end

# ╔═╡ ad2b0861-50db-487f-b0fd-b6f95fa0358f
md"""
## 2. Phase Distortions in Scattered Wavefields
 
When acoustic or elastic waves propagate through media with small-scale heterogeneities, such as rough interfaces or sub-wavelength inclusions, each scattering event imposes a local phase distortion. The cumulative effect is a frequency-dependent, spatially varying phase perturbation that differs from record to record, decorrelating otherwise coherent recordings.
 
Scattering-induced phase randomization is the physical origin of speckle⁸, observed across acoustics, radar, ultrasound, and interferometric remote sensing. In the complex signal domain, speckle manifests as multiplicative noise: the observed field is the product of a coherent signal and a random phasor whose phase approaches a uniform distribution under strong scattering. Analogous multiplicative phase distortions have recently been identified in broadband seismic data, where near-surface scattering introduces random, frequency-dependent phase perturbations across recorded signals, manifesting as severe wavelet distortions in the time domain⁹.
 
Conceptually, a clean spectrum
"""

# ╔═╡ 3d3e9271-505e-42e9-9fac-52b54a02c8f4
md"""
```math
A_{\text{clean}}(\omega) = |A_{\text{clean}}(\omega)|\, e^{i\theta_{\text{clean}}(\omega)} \qquad (1)
```
"""

# ╔═╡ 5073c1d5-9b7c-4636-8fa1-3c8882ac30b6
md"""
acquires a record-dependent phase perturbation ``\delta_k(\omega)``. The perturbations ``\delta_k(\omega)`` are statistically independent across records and frequencies, and frequency-dependent in character. The spectral amplitude ``|A_{\text{clean}}(\omega)|`` remains unchanged; the decorrelation is entirely phase driven.
 
The perturbed spectrum of the ``k``-th recording in a multichannel seismic ensemble can therefore be written as
"""

# ╔═╡ 19ea21f9-32b3-4e10-9439-d1731a0f5681
md"""
```math
A_k(\omega) = |A_{\text{clean}}(\omega)|\, e^{i\theta_k(\omega)}, \qquad \theta_k(\omega) = \theta_{\text{clean}}(\omega) + \delta_k(\omega). \qquad (2)
```
"""

# ╔═╡ 5b2901db-9918-4bda-a084-0e4ada7e0301
md"""
Here, ``k`` indexes individual recorded traces in the ensemble. In the time domain, this phase perturbation produces distorted waveform shapes (Fig. 1c), mimicking scattering-induced phase scrambling observed in field data.
 
In practice, trace amplitudes also vary due to geometric spreading, attenuation, source coupling, and scattering losses. Incorporating trace-dependent amplitude scaling ``\alpha_k`` drawn from a uniform distribution on ``[0, 1]`` gives
"""

# ╔═╡ 28465d38-af5e-4b70-b6b2-85785fc7c1cb
md"""
```math
A'_k(\omega) = \alpha_k\, |A_{\text{clean}}(\omega)|\, e^{i\theta_k(\omega)}. \qquad (3)
```
"""

# ╔═╡ 74b5f360-988d-4f79-a96c-83682a02d2ae
md"""
This two-stage synthetic design: phase-only perturbations followed by combined phase and amplitude variability, enables controlled evaluation under idealized (uniform-amplitude) and realistic (variable-amplitude) conditions (Fig. 1d).
 
Suppressing these perturbations requires an operator that acts exclusively on the phase component of the spectrum.
"""

# ╔═╡ 40c12c7b-18be-4cf5-a42b-0ee57f30281b
t = range(start=0, stop=length(klauder)*T, length=length(klauder))

# ╔═╡ a2d7634f-9bcb-44cd-9f98-cb582cad0e3f
function trace_spectrum(data, dt; extras=false)
	nt, ntr = size(data)
	freq = fftfreq(nt, 1/dt)
	pos_inds = findall(x -> x > 0, freq) # excludes DC
	nfreq = length(pos_inds)
	phase_vals = zeros(Float32, nfreq, ntr)
	amplitude_vals = zeros(Float32, nfreq, ntr)
	f_zero = zeros(ComplexF64, 1, ntr)
	for k in 1:ntr
		fft_trace = fft(data[:, k])
		f_zero[1, k] = fft_trace[1]
		phase_vals[:, k] = angle.(fft_trace[pos_inds])
		amplitude_vals[:, k] = abs.(fft_trace[pos_inds])
	end
	return extras ? (phase_vals, amplitude_vals, freq[pos_inds], f_zero) :
	(phase_vals, amplitude_vals)
end

# ╔═╡ fb919903-e6c4-45a4-97bd-b7e64b23c1ac
phase_clean, amplitude_clean, positive_freq, f_zero =
trace_spectrum(clean_traces, T, extras=true)

# ╔═╡ 9324270e-ca58-46d9-a2ad-bbb7d259cf26
function perturb_phase(phase_clean, amplitude_clean, f_zero, num_wiggles;
                        kappa_start=0.9, kappa_stop=0.1, seed=2025)

    Random.seed!(seed)

    n_freqs      = size(phase_clean, 1)
    kappa_values = range(start=kappa_start, stop=kappa_stop, length=n_freqs)
    var_values   = [1 - besseli(1, k)/besseli(0, k) for k in kappa_values]

    # apply phase perturbations
    phase_perturbed = similar(phase_clean)
    for i in 1:n_freqs
        pp             = phase_clean[i, :]
        kappat         = kappa_values[i]
        von_mises_dist = Distributions.VonMises(0.0, kappat)
        random         = rand(von_mises_dist, size(phase_clean, 2))
        phase_perturbed[i, :] .= pp .+ random
    end
    phase_perturbed = mod.(phase_perturbed .+ π, 2π) .- π

    # reconstruct perturbed signal via inverse FFT
    new_f = amplitude_clean .* exp.(1.0im .* phase_perturbed)
    perturbed_signal_half = zeros(ComplexF64, size(new_f, 1), num_wiggles)
    for i in 1:size(new_f, 2)
        perturbed_signal_half[:, i] = reverse(conj(new_f[:, i]))
    end

    perturbed_signal = vcat(f_zero, new_f, perturbed_signal_half)
    ifft_perturbed   = zeros(Float64, size(perturbed_signal, 1), num_wiggles)
    for i in 1:size(perturbed_signal, 2)
        ifft_perturbed[:, i] = real(ifft(perturbed_signal[:, i]))
    end

    return ifft_perturbed, phase_perturbed, collect(var_values)
end

# ╔═╡ 3519b3d0-fa46-494b-9f38-70d6aeae0b82
phase_perturbed_traces, phase_perturbed, var_imposed = perturb_phase(
        phase_clean, amplitude_clean, f_zero, num_wiggles;
        kappa_start=0.9, kappa_stop=0.1, seed=2025)

# ╔═╡ 0173fbd2-044d-4eff-bbc0-e3c3ec67fa88
function plot_wiggles_va(traces, t;
                          n_traces=20, step=1, scale=1.0, lw=1.8,
                          title="", xlabel="Trace Index", ylabel="Time (s)",
                          xtick_labels=["0", "N"])

    traces_subset = traces[:, 1:min(n_traces, size(traces, 2))]
    cols  = 1:step:size(traces_subset, 2)
    tseg  = collect(t)
    nplot = length(cols)

    p = plot(legend=false, grid=false, yflip=true,
             title=title, xlabel=xlabel, ylabel=ylabel,
             size=(1000, 300),
             fontsize=12, xtickfontsize=12, ytickfontsize=12,
             guidefontsize=12, titlefontsize=12, legendfontsize=12)

    for (j, k) in enumerate(cols)
        trc  = traces_subset[:, k]
        xwig = trc .* scale .+ j

        plot!(p, xwig, tseg, color=:black, lw=lw)

        mask   = trc .> 0
        edges  = findall(diff([false; mask; false]) .!= 0)
        starts = edges[1:2:end]
        stops  = edges[2:2:end] .- 1

        for (s, e) in zip(starts, stops)
            xpoly = vcat(xwig[s:e], fill(j, e - s + 1))
            ypoly = vcat(tseg[s:e], reverse(tseg[s:e]))
            plot!(p, xpoly, ypoly, seriestype=:shape,
                  fillcolor=:black, linecolor=:black, lw=0)
        end
    end

    xlims!(p, 0.5, nplot + 0.5)
    xticks!(p, [1, nplot], xtick_labels)

    return p
end

# ╔═╡ 44737414-985a-4cf5-80e5-ba822bb16231
function phase_masking_via_local_stack(perturbed_signal, amplitude, f_zero, fs;
                                        local_window=6)

    n_samples, n_traces = size(perturbed_signal)
    n_windows = n_traces - local_window + 1

    # ── local stacking 
    local_stack = Matrix{Float64}(undef, n_samples, n_windows)
    for j in 1:n_windows
        col_indices = j:(j + local_window - 1)
        local_stack[:, j] = vec(sum(perturbed_signal[:, col_indices], dims=2) / local_window)
    end

    # ── extract predicted phases from local stack 
    predicted_phases = zeros(Float64, 86, n_windows)
    for i in 1:n_windows
        f       = fft(local_stack[:, i])
        freq    = fftfreq(length(f), fs)
        pos_idx = findall(x -> x > 0, freq)
        predicted_phases[:, i] = angle.(f[pos_idx])
    end

    # ── reconstruct phase-corrected signal 
    new_f = amplitude[:, 1:n_windows] .* exp.(1.0im .* predicted_phases)
    perturbed_signal_half = zeros(ComplexF64, size(new_f, 1), n_windows)
    for i in 1:n_windows
        perturbed_signal_half[:, i] = reverse(conj(new_f[:, i]))
    end

    perturbed_signal_full = vcat(f_zero[1:n_windows]', new_f, perturbed_signal_half)
    ifft_output = zeros(Float64, size(perturbed_signal_full, 1), n_windows)
    for i in 1:n_windows
        ifft_output[:, i] = real(ifft(perturbed_signal_full[:, i]))
    end

    return ifft_output, predicted_phases, local_stack
end

# ╔═╡ 85e9ae48-a909-4c05-b8b5-d8d1ce39329f
pls_traces, predicted_phases_ls, local_stack = phase_masking_via_local_stack(
        phase_perturbed_traces, amplitude_clean, f_zero, fs; local_window=6)

# ╔═╡ f0b9b033-5893-4e2e-8c72-d94547f50b11
function circular_mean(α; w = ones(size(α)), dims = 1)
r = sum(w .* cis.(α); dims)
μ = angle.(r)
μ = length(μ) == 1 ? μ[1] : μ
return μ
end

# ╔═╡ 2b0bebc0-7159-4446-baec-ac5a5e5bb784
function phase_stats(phase_matrix)
	nfreq = size(phase_matrix, 1)
	mean_vals = zeros(Float64, nfreq)
	var_vals = zeros(Float64, nfreq)
	for i in 1:nfreq
	row = phase_matrix[i, :]
	mean_vals[i] = circular_mean(row)
	var_vals[i] = 1 - circ_r(row)
	end
	return mean_vals, var_vals
end

# ╔═╡ 3932ee58-1132-465a-9b03-33e690cecb81
function phase_masking_via_circular_mean(perturbed_phases, amplitude, f_zero, num_wiggles;
                                          local_window=6)

    half_window = div(local_window, 2)
    n_freqs, n_traces = size(perturbed_phases)

    # ── compute circular mean phases 
    predicted_phases = Matrix{Float64}(undef, n_freqs, n_traces)
    for j in 1:n_traces
        start_idx = max(1, j - half_window)
        end_idx   = min(n_traces, j + half_window - 1)
        mid_idx   = min(n_traces, max(1, div(start_idx + end_idx, 2)))

        window_phases    = perturbed_phases[:, start_idx:end_idx]
        mean_vals, _     = phase_stats(window_phases)

        predicted_phases[:, j]       .= mean_vals
        predicted_phases[:, mid_idx] .= mean_vals
    end

    # ── reconstruct phase-corrected signal 
    new_f = amplitude .* exp.(1.0im .* predicted_phases)
    perturbed_signal_half = zeros(ComplexF64, size(new_f, 1), num_wiggles)
    for i in 1:size(new_f, 2)
        perturbed_signal_half[:, i] = reverse(conj(new_f[:, i]))
    end

    perturbed_signal_full = vcat(f_zero, new_f, perturbed_signal_half)
    ifft_output = zeros(Float64, size(perturbed_signal_full, 1), num_wiggles)
    for i in 1:size(perturbed_signal_full, 2)
        ifft_output[:, i] = real(ifft(perturbed_signal_full[:, i]))
    end

    return ifft_output, predicted_phases
end

# ╔═╡ cebefee6-4653-48d9-893d-b221b0e83535
pcm_traces, predicted_phases_pcm = phase_masking_via_circular_mean(
        phase_perturbed, amplitude_clean, f_zero, num_wiggles; local_window=6)


# ╔═╡ 31479529-c0cf-42f8-8fea-465084fe8c10
begin
    f   = 2
    idx = 1:4995

    # 1a — clean traces
    p1a = plot_wiggles_va(clean_traces, t;
                          n_traces=50, step=1, scale=1.5, lw=1.0,
                          title="(a) Clean seismic records",
                          xtick_labels=["0", "5000"])

    # 1c — uniform amplitudes perturbed phases
    p1c = plot_wiggles_va(phase_perturbed_traces, t;
                          n_traces=100, step=1, scale=1.5, lw=1.0,
                          title="(c) Uniform amplitudes and perturbed phases",
                          xtick_labels=["0", "5000"])

    # 1e — PLS uniform amplitudes
    p1e = plot_wiggles_va(pls_traces, t;
                          n_traces=100, step=1, scale=1.5, lw=1.0,
                          title="(e) Phase masking via local stack (PLS)",
                          xtick_labels=["0", "5000"])

    # 1g — PCM uniform amplitudes
    p1g = plot_wiggles_va(pcm_traces, t;
                          n_traces=100, step=1, scale=1.5, lw=1.0,
                          title="(g) Phase masking via circular mean (PCM)",
                          xtick_labels=["0", "5000"])

    # 1i — phase comparison uniform amplitudes
    p1i = plot(phase_perturbed[f, idx],
               label="Perturbed", color=:red, linestyle=:solid, lw=2,
               size=(1000, 300), grid=false,
               ylabel="Phase", xlabel="Trace Index",
               title="(i)",
               fontsize=12, labelfontsize=12, legendfontsize=12,
               tickfontsize=12, xguidefontsize=12, yguidefontsize=12)
    plot!(p1i, predicted_phases_ls[f, idx],  label="PLS",   color=:blue,     lw=2)
    plot!(p1i, predicted_phases_pcm[f, idx], label="PCM",   color="#6fa76f", lw=2)
    plot!(p1i, phase_clean[f, idx],          label="Clean", color=:black,    lw=2)

    plot(p1a, p1c, p1e, p1g, p1i,
         layout=(5, 1),
         size=(1000, 1500),
         left_margin=8Plots.mm,
         bottom_margin=6Plots.mm,
         top_margin=4Plots.mm,
         dpi=300)
end

# ╔═╡ 07598432-1156-4d27-89c6-422dbaef6914
function scale_amplitudes(amplitude_clean, num_wiggles;
                           scale_min=0.0, scale_max=1.0, seed=2025)
    Random.seed!(seed)
    scales     = scale_min .+ (scale_max - scale_min) .* rand(num_wiggles)
    amp_scaled = amplitude_clean .* scales'
    return amp_scaled, scales
end

# ╔═╡ ed98073b-8ac5-4bea-b0d6-cf646039196c
amp_scaled, scales = scale_amplitudes(amplitude_clean, num_wiggles;
                                          scale_min=0.0, scale_max=1.0, seed=2025)

# ╔═╡ 5d44094f-b541-45f3-b164-a43af866d911
scaled_perturbed_traces, phase_perturbed_scaled, var_imposed_scaled = perturb_phase(
        phase_clean, amp_scaled, f_zero, num_wiggles;
        kappa_start=0.9, kappa_stop=0.1, seed=2025)

# ╔═╡ dac97421-f42b-41e7-8569-a2cb5c5f2085
pls_traces_scaled, predicted_phases_ls_scaled, local_stack_scaled =
        phase_masking_via_local_stack(scaled_perturbed_traces, amp_scaled, f_zero, fs; local_window=6)

# ╔═╡ f8b66581-f13a-451d-a77a-c19395d89ff1
pcm_traces_scaled, predicted_phases_pcm_scaled =
        phase_masking_via_circular_mean(phase_perturbed_scaled, amp_scaled, f_zero, num_wiggles;
                                         local_window=6)

# ╔═╡ b3f6f24f-0260-4d0b-b80b-e60d5ffd2b3a
begin
    # 1b — amplitude comparison
    p1b = plot(amplitude_clean[2, 1:5000], label="Clean", color=:black, lw=2,
               title="(b) Amplitude scale factors",
               xlabel="Trace Index", ylabel="Amplitude",
               grid=false, legend=false,
               fontsize=12, labelfontsize=12, legendfontsize=12,
               tickfontsize=12, xguidefontsize=12, yguidefontsize=12,
               size=(1000, 300), dpi=300)
    plot!(p1b, amp_scaled[2, 1:5000], label="Scaled",
          color=:orange, linestyle=:dash, lw=0.7)

    # 1d — variable amplitudes with perturbed phases
    p1d = plot_wiggles_va(scaled_perturbed_traces, t;
                          n_traces=100, step=1, scale=1.5, lw=1.0,
                          title="(d) Variable amplitudes with perturbed phases",
                          xtick_labels=["0", "5000"])

    # 1f — PLS on scaled amplitudes
    p1f = plot_wiggles_va(pls_traces_scaled, t;
                          n_traces=100, step=1, scale=1.5, lw=1.0,
                          title="(f) Phase masking via local stack (PLS)",
                          xtick_labels=["0", "5000"])

    # 1h — PCM on scaled amplitudes
    p1h = plot_wiggles_va(pcm_traces_scaled, t;
                          n_traces=100, step=1, scale=1.5, lw=1.0,
                          title="(h) Phase masking via circular mean (PCM)",
                          xtick_labels=["0", "5000"])

    # 1j — phase comparison scaled amplitudes
    p1j = plot(phase_perturbed_scaled[f, idx],
               label="Perturbed", color=:red, linestyle=:solid, lw=2,
               size=(1000, 300), grid=false,
               ylabel="Phase", xlabel="Trace Index",
               title="(j)",
               fontsize=12, labelfontsize=12, legendfontsize=12,
               tickfontsize=12, xguidefontsize=12, yguidefontsize=12)
    plot!(p1j, predicted_phases_ls_scaled[f, idx],  label="PLS",   color=:blue,     lw=2)
    plot!(p1j, predicted_phases_pcm_scaled[f, idx], label="PCM",   color="#6fa76f", lw=2)
    plot!(p1j, phase_clean[f, idx],                 label="Clean", color=:black,    lw=2)

    plot(p1b, p1d, p1f, p1h, p1j,
         layout=(5, 1),
         size=(1000, 1500),
         left_margin=8Plots.mm,
         bottom_margin=6Plots.mm,
         top_margin=4Plots.mm,
         dpi=300)
end

# ╔═╡ 3d6b96fd-53e9-4d89-8327-de4429cb83d8
md"""
## 3. Phase-Only Coherence Conditioning via Phase Masking
 
Restoring coherence across a scattered wavefield requires correcting record-dependent phase perturbations. We focus exclusively on phase distortions, which are often overlooked in conventional seismic signal-processing workflows. Amplitude effects are not the scope of this letter; instead, our goal is to correct phase inconsistencies while preserving the original spectral amplitudes.
 
A phase mask accomplishes this as a multiplicative, unit-magnitude operator in the complex frequency domain that replaces the phase of each record with a locally estimated representative phase while leaving spectral amplitudes exactly unchanged.
 
For record ``k``, the phase mask at frequency ``\omega`` is defined as
"""

# ╔═╡ e91852af-5b48-41ca-b3ff-51d127d3d7e1
md"""
```math
PM_k(\omega) = \exp\!\left(i\left[\hat{\theta}_k(\omega) - \theta_k(\omega)\right]\right). \qquad (4)
```
"""

# ╔═╡ c623aa1c-ebcf-4da4-a2c6-f417fb6abb75
md"""
The mask is applied directly to the original complex spectrum and the masked output is
"""

# ╔═╡ 65e06a4c-287c-4a17-822b-a9d71de61368
md"""
```math
\tilde{A}_k(\omega) = |A_k(\omega)|\, e^{i\hat{\theta}_k(\omega)}. \qquad (5)
```
"""

# ╔═╡ 4b5f5011-156a-421e-ba81-e526f3d2a79e
md"""
where ``\hat{\theta}_k(\omega)`` is the representative phase estimated from local records that sample comparable propagation paths. By construction ``|PM_k(\omega)| = 1`` at all frequencies, so spectral energy is preserved and all changes solely arise from phase conditioning.
 
The two methods compared next share this formulation and differ only in how ``\hat{\theta}_k(\omega)`` is estimated from the local ensemble.
 
### Phase Masking via Local Stack (PLS)
 
One approach to estimating the representative phase is to average the ensemble of ``n`` recorded seismic signals ``x_k(t)`` in the time domain. The phase of the resulting stack,
"""

# ╔═╡ d4ae36d7-efa8-497e-b581-5268c0e3f863
md"""
```math
x_{\text{stack}}(t) = \frac{1}{n}\sum_{k=1}^{n} x_k(t), \qquad (6)
```
"""

# ╔═╡ 10eab8aa-247e-424d-81d3-e6feb0832d47
md"""
defines ``\theta_{\text{stack}}(\omega)`` as the phase spectrum of ``x_{\text{stack}}(t)``, which then serves as the representative phase: ``\hat{\theta}_{PLS}(\omega) = \theta_{\text{stack}}(\omega)``. PLS is reliable when trace amplitudes are uniform. However, the time-domain stack is amplitude weighted, so ``\hat{\theta}_{PLS}(\omega)`` becomes biased whenever amplitude variability is present.
 """

# ╔═╡ 8258d07d-4a86-4965-a509-e1645ca25ced
md"""
 
### Phase Masking via Circular Mean (PCM)
 
PCM addresses the amplitude bias inherent in PLS by removing amplitude information before the ensemble average is formed. Each spectrum is normalized to unit magnitude, ensuring all traces contribute equally regardless of their spectral energy. The representative phase for frequency ``\omega`` and ensemble of ``n`` records is estimated using the circular mean¹⁰,
"""

# ╔═╡ f05c7e7f-4d99-444f-b179-d79375392896
md"""
```math
\hat{\theta}_{PCM}(\omega) = \arg\!\left(\frac{1}{n}\sum_{k=1}^{n}\frac{A_k(\omega)}{|A_k(\omega)|}\right). \qquad (7)
```
"""

# ╔═╡ 2cca3d97-5696-4fbe-ad34-368cbb4aa0b4
md"""
Since normalization precedes averaging, all traces contribute equally regardless of their spectral energy. ``\hat{\theta}_{PCM}(\omega)`` is amplitude-invariant by construction. The circular mean provides a robust estimate of the dominant signal phase direction¹¹.
 
The two estimators are equivalent under uniform amplitudes and systematically diverge as amplitude variability increases. PLS is therefore a limiting case of the general phase-masking framework: its amplitude weighting is harmless when trace energies are approximately equal but introduces a bias toward high-amplitude records otherwise. PCM eliminates this bias by construction, ensuring equal contribution from all traces regardless of spectral energy.
 
The practical consequence is visible in Figs. 1e–1h, which show the conditioned records under both idealized (uniform-amplitude) and realistic (variable-amplitude) conditions. Under uniform amplitudes, PLS and PCM produce indistinguishable results (Figs. 1e, 1g), consistent with the theoretical equivalence in that limit. Under variable amplitudes, PCM yields noticeably more consistent phase alignment across the ensemble (Fig. 1h), while the PLS output retains residual phase scatter attributable to amplitude bias (Fig. 1f).
"""

# ╔═╡ 65b1a5ba-c4e2-493d-a5cf-ec595fba83d7
wrap_phase_local(x) = mod.(x .+ π, 2π) .- π

# ╔═╡ 4a52c424-bc43-436c-9fc6-0bd4bb978a2a
function make_phasor_circle(θ_pert, θ_ls, θ_pc, θ_clean; title_str="")
        n = length(θ_clean)
        r = range(0.05, 1.0; length=n)
        rmax = maximum(r)

        x_clean = cos.(θ_clean) .* r;  y_clean = sin.(θ_clean) .* r
        x_pert  = cos.(θ_pert)  .* r;  y_pert  = sin.(θ_pert)  .* r
        x_pc    = cos.(θ_pc)    .* r;  y_pc    = sin.(θ_pc)    .* r
        x_ls    = cos.(θ_ls)    .* r;  y_ls    = sin.(θ_ls)    .* r

        unit_x = cos.(range(-π, π, length=400))
        unit_y = sin.(range(-π, π, length=400))

        p = plot(unit_x .* rmax, unit_y .* rmax,
                 color=:black, lw=2, label="",
                 aspect_ratio=:equal,
                 axis=false, ticks=false, grid=false,
                 framestyle=:none, size=(600, 600),
                 title=title_str)

        plot!(p, x_pert,  y_pert,  color=:red,      lw=1, ls=:solid, label="Perturbed")
        plot!(p, x_ls,    y_ls,    color=:blue,     lw=1, ls=:solid, label="PLS")
        plot!(p, x_pc,    y_pc,    color="#6fa76f", lw=1, ls=:solid, label="PCM")
        plot!(p, x_clean, y_clean, color=:black,    lw=2,            label="Clean")

        # annotations
        annotate!(p, -rmax * 1.05, 0.0, text("-π,\nπ", :black, 8, :right))
        annotate!(p,  rmax * 1.05, 0.0, text("0",      :black, 8, :left))

        return p
end

# ╔═╡ 4ea2887e-b302-4377-bc4f-b838072576e4
begin
    number_of_segments(θ₁, θ₂) = max(2, round(Int, abs(θ₂ - θ₁) * 90 / π))

    arc(θ₁, θ₂, r=1, x=0, y=0) = Tuple{Float64,Float64}[
        (x + r*cos(θ), y + r*sin(θ))
        for θ in range(θ₁, stop=θ₂, length=number_of_segments(θ₁, θ₂))
    ]

    circle(r, x=0, y=0) = Shape(arc(0, 2π, r, x, y))

    sector(θ₁, θ₂, r=1, x=0, y=0)::Shape =
        Shape(vcat(arc(θ₁, θ₂, r, x, y), Tuple{Float64,Float64}((x, y))))

    function fit_hist(a, binwidth, degrees::Bool, axial::Bool;
                      weights=ones(eltype(a), length(a)))
        binwidth isa Bool && throw(ArgumentError("`binwidth` supplied as a Bool"))
        length(weights) == length(a) ||
            throw(ArgumentError("`weights` must have same length as number of angles"))
        n    = degrees ? round(Int, 360/binwidth) : round(Int, 2π/binwidth)
        n    = max(1, n)
        bins = range(0, stop=2π, length=n+1)
        data = degrees ? deg2rad.(mod.(a, 360)) : mod.(a, 2π)
        data .+= 10*eps(float(eltype(a)))
        axial && (data .= mod.(data, π))
        h = StatsBase.fit(StatsBase.Histogram, data,
                          StatsBase.FrequencyWeights(weights), bins, closed=:left)
        axial && (h.weights[end÷2+1:end] .= h.weights[1:end÷2])
        return h
    end
end

# ╔═╡ 3e18f298-6ba0-472f-b421-c1bf043b5080
function cplot_histogram_overlay3(a1, a2, a3, binwidth;
    degrees::Bool=false, azimuth::Bool=false, axial::Bool=false,
    weights1=ones(eltype(a1), length(a1)),
    weights2=ones(eltype(a2), length(a2)),
    weights3=ones(eltype(a3), length(a3)),
    normalize::Symbol=:none,
    circ=(lc=:black, fill=nothing),
    style1=(; fillcolor=:white, linecolor="#6fa76f", fillalpha=0.25, linewidth=1),
    style2=(; fillcolor=:white, linecolor=:blue,     fillalpha=0.25, linewidth=1),
    style3=(; fillcolor=:white, linecolor=:red,      fillalpha=0.25, linewidth=1),
    title_str=""
)
    h1 = fit_hist(a1, binwidth, degrees, axial; weights=weights1)
    h2 = fit_hist(a2, binwidth, degrees, axial; weights=weights2)
    h3 = fit_hist(a3, binwidth, degrees, axial; weights=weights3)

    edges = collect(h1.edges[1])
    n = length(h1.weights)
    (length(h2.weights) == n && length(h3.weights) == n) ||
        throw(ArgumentError("Histograms ended up with different bin counts; use the same binwidth."))

    edges_use = azimuth ? (π/2 .- edges) : edges

    r1 = Float64.(h1.weights)
    r2 = Float64.(h2.weights)
    r3 = Float64.(h3.weights)

    if normalize != :none
        s1, s2, s3 = sum(r1), sum(r2), sum(r3)
        (s1 > 0 && s2 > 0 && s3 > 0) || throw(ArgumentError("One histogram has zero total weight."))
        r1 ./= s1; r2 ./= s2; r3 ./= s3
        if normalize == :density
            bw = edges[2] - edges[1]
            r1 ./= bw; r2 ./= bw; r3 ./= bw
        elseif normalize != :prob
            throw(ArgumentError("normalize must be :none, :prob, or :density"))
        end
    end

    R = max(maximum(r1), maximum(r2), maximum(r3))
    R == 0 && (R = 1.0)

    p = plot(aspect_ratio=:equal, xlim=(-R, R), ylim=(-R, R),
             legend=false, showaxis=false, grid=false, title=title_str)

    circ !== nothing && plot!(p, circle(R); circ...)

    plot!(p, sector.(edges_use[1:n], edges_use[2:n+1], r1);
          fill=style1.fillcolor, fillalpha=style1.fillalpha,
          lc=style1.linecolor, lw=style1.linewidth)
    plot!(p, sector.(edges_use[1:n], edges_use[2:n+1], r2);
          fill=style2.fillcolor, fillalpha=style2.fillalpha,
          lc=style2.linecolor, lw=style2.linewidth)
    plot!(p, sector.(edges_use[1:n], edges_use[2:n+1], r3);
          fill=style3.fillcolor, fillalpha=style3.fillalpha,
          lc=style3.linecolor, lw=style3.linewidth)

    return p
end

# ╔═╡ 46a4c8ea-33d0-46a8-aa63-4c3f1c04fad9
md"""
## 4. Quantitative Assessment of Phase Conditioning: PLS vs PCM
 
Having introduced both estimators, we now evaluate their performance under synthetic conditions using two measures: circular variance and SNR.
 
The degree of phase coherence is quantified using circular variance,
"""

# ╔═╡ 4a115cf9-648f-4719-8d87-7d8e4796bdeb
md"""
```math
V(\omega) = 1 - \left|\frac{1}{n}\sum_{k=1}^{n}\frac{A_k(\omega)}{|A_k(\omega)|}\right|, \qquad (8)
```
"""

# ╔═╡ 9bcaedf0-c58d-439e-933c-4761439bd1a4
md"""
which ranges from 0 (perfectly coherent) to 1 (perfectly uniform), providing an objective, amplitude-independent measure of phase dispersion across the ensemble.
 
Signal-to-noise ratio (SNR) is estimated using semblance ``S`` as an intermediate variable¹²,
"""

# ╔═╡ c755f57d-0b71-4e06-81ca-7b64e581d544
md"""
```math
\text{SNR}_{dB} = 10\log_{10}\frac{S}{1-S}. \qquad (9)
```
"""

# ╔═╡ 36ea57c3-1ad5-473d-b2c1-15431754b670
begin
    binw   = 2π/70
	
    # uniform amplitude stats
    var_pert_u = round(1 - circ_r(phase_perturbed[f, idx]),      digits=3)
    var_ls_u   = round(1 - circ_r(predicted_phases_ls[f, idx]),  digits=3)
    var_pc_u   = round(1 - circ_r(predicted_phases_pcm[f, idx]), digits=3)

    # variable amplitude stats
    var_pert_s = round(1 - circ_r(phase_perturbed_scaled[f, idx]),      digits=3)
    var_ls_s   = round(1 - circ_r(predicted_phases_ls_scaled[f, idx]),  digits=3)
    var_pc_s   = round(1 - circ_r(predicted_phases_pcm_scaled[f, idx]), digits=3)
end

# ╔═╡ 0573f916-bbb6-40e6-a6bd-8aea25ff36f1
md"""
We compare PLS and PCM under two scenarios designed to isolate the conditions under which the methods agree and where they differ: uniform amplitudes across records and variable amplitudes across records.
"""

# ╔═╡ f8341612-b575-4e9c-bbdb-424c405664c8
md""" 
### Case 1: Phase Perturbations Only with Uniform Amplitudes
 
Under uniform amplitudes, both PLS and PCM reduce phase variability equivalently, a result that confirms that the two estimators converge in the idealized limit. At 5 Hz, the perturbed ensemble of 5000 records exhibits a circular variance of ``V = 0.597`` (Fig. 2c), indicating broad phase scatter around the dominant phase direction (Fig. 2a). After applying either PLS or PCM, variance drops identically to ``V = 0.268`` (Fig. 2d–e), a 55% reduction, reflecting tighter concentration of the phase distribution around the circular mean (Fig. 2b, blue and green curves overlap). Consistent with this, SNR improves from ``-8.09`` dB to ``-0.55`` dB after PLS and ``-0.56`` dB after PCM. This equivalence holds across the full bandwidth (Fig. 3a).
 
When amplitudes are uniform, the time-domain stack provides an unbiased estimate of signal phase, making PLS a computationally simpler choice.
"""

# ╔═╡ a4b60ad9-dd61-4980-bab7-373d7d6ba466
begin
    θ_clean_u = wrap_phase_local(phase_clean[f, idx])
    θ_pert_u  = wrap_phase_local(phase_perturbed[f, idx])
    θ_ls_u    = wrap_phase_local(predicted_phases_ls[f, idx])
    θ_pc_u    = wrap_phase_local(predicted_phases_pcm[f, idx])

    p2a = make_phasor_circle(θ_pert_u, θ_ls_u, θ_pc_u, θ_clean_u;
                              title_str="(a) Uniform amplitudes (5 Hz)")

    p2b = cplot_histogram_overlay3(
        predicted_phases_pcm[f, idx],
        predicted_phases_ls[f, idx],
        phase_perturbed[f, idx],
        binw;
        normalize=:none,
        circ=(lc=:black, fill=nothing),
        title_str="(b) Uniform\nPCM(V=$var_pc_u), PLS(V=$var_ls_u), Pert(V=$var_pert_u)")

    plot(p2a, p2b, layout=(1,2), size=(1200,500),
         left_margin=8Plots.mm, bottom_margin=6Plots.mm, dpi=300)
end

# ╔═╡ b3bfcd66-2959-4192-a44b-4f27583ffe0d
begin
    p2c = histogram(phase_perturbed[f, idx], bins=70, normalize=true,
                    fill=:white, lc=:red, legend=false, grid=false,
                    title="(c) Perturbed\n(V=$var_pert_u)", showaxis=false)
    ylims!(p2c, 0, 0.7)
    p2d = histogram(predicted_phases_ls[f, idx], bins=70, normalize=true,
                    fill=:white, lc=:blue, legend=false, grid=false,
                    title="(d) PLS\n(V=$var_ls_u)", showaxis=false)
    ylims!(p2d, 0, 0.7)
    p2e = histogram(predicted_phases_pcm[f, idx], bins=70, normalize=true,
                    fill=:white, lc="#6fa76f", legend=false, grid=false,
                    title="(e) PCM\n(V=$var_pc_u)", showaxis=false)
    ylims!(p2e, 0, 0.7)
    plot(p2c, p2d, p2e, layout=(1,3), size=(1200,400),
         left_margin=6Plots.mm, bottom_margin=6Plots.mm, dpi=300)
end

# ╔═╡ cdf2bd79-8064-488a-8eb5-d4690fc4e560
md"""
 
### Case 2: Phase Perturbations with Variable Amplitudes
 
Introducing record-dependent amplitude variability breaks the equivalence: PLS phase estimates become biased by amplitude mixing, whereas PCM remains robust by construction. This divergence is both statistically significant and physically consequential. At 5 Hz, after applying PLS, variance decreases only to ``V = 0.32`` (Fig. 2i), a partial reduction of 45%, because the time-domain stack is dominated by high-amplitude records, pulling the representative phase away from the true signal phase. PCM eliminates this bias and reduces variance to ``V = 0.26`` (Fig. 2j). The circular representations in Figs. 2f–g make this contrast visually explicit: after PCM, phase clusters are as tight as in Case 1; after PLS, they remain substantially broader. Fig. 3b shows the full spectrum: PCM maintains consistently lower phase variance than PLS. This is also reflected in SNR: starting from the same baseline of ``-9.14`` dB, PLS achieves ``-3.60`` dB (a gain of 5.54 dB), while PCM reaches ``-2.52`` dB (a gain of 6.62 dB).
 
These results are summarized in Table 1. Together, they establish that PLS and PCM are interchangeable only in the idealized uniform-amplitude regime. Once amplitude variability is present — the realistic condition in any field acquisition, explicit amplitude-phase decoupling via circular-mean phase masking is essential for reliable coherence conditioning.
"""

# ╔═╡ 07322a04-dfeb-4c12-97e8-b3e993baf352
# Cell for p2f and p2g
begin
    θ_clean_s = wrap_phase_local(phase_clean[f, idx])
    θ_pert_s  = wrap_phase_local(phase_perturbed_scaled[f, idx])
    θ_ls_s    = wrap_phase_local(predicted_phases_ls_scaled[f, idx])
    θ_pc_s    = wrap_phase_local(predicted_phases_pcm_scaled[f, idx])

    p2f = make_phasor_circle(θ_pert_s, θ_ls_s, θ_pc_s, θ_clean_s;
                              title_str="(f) Variable amplitudes (5 Hz)")

    p2g = cplot_histogram_overlay3(
        predicted_phases_pcm_scaled[f, idx],
        predicted_phases_ls_scaled[f, idx],
        phase_perturbed_scaled[f, idx],
        binw;
        normalize=:none,
        circ=(lc=:black, fill=nothing),
        title_str="(g) Variable\nPCM(V=$var_pc_s), PLS(V=$var_ls_s), Pert(V=$var_pert_s)")

    plot(p2f, p2g, layout=(1,2), size=(1200,500),
         left_margin=8Plots.mm, bottom_margin=6Plots.mm, dpi=300)
end

# ╔═╡ 1948cc93-38ca-4be4-893a-531827b27312
begin
    p2h = histogram(phase_perturbed_scaled[f, idx], bins=70, normalize=true,
                    fill=:white, lc=:red, legend=false, grid=false,
                    title="(h) Perturbed\n(V=$var_pert_s)", showaxis=false)
    ylims!(p2h, 0, 0.7)
    p2i = histogram(predicted_phases_ls_scaled[f, idx], bins=70, normalize=true,
                    fill=:white, lc=:blue, legend=false, grid=false,
                    title="(i) PLS\n(V=$var_ls_s)", showaxis=false)
    ylims!(p2i, 0, 0.7)
    p2j = histogram(predicted_phases_pcm_scaled[f, idx], bins=70, normalize=true,
                    fill=:white, lc="#6fa76f", legend=false, grid=false,
                    title="(j) PCM\n(V=$var_pc_s)", showaxis=false)
    ylims!(p2j, 0, 0.7)
    plot(p2h, p2i, p2j, layout=(1,3), size=(1200,400),
         left_margin=6Plots.mm, bottom_margin=6Plots.mm, dpi=300)
end

# ╔═╡ a506d20f-6b24-476c-9a46-8301c75c4b05
md"""
**Table 1.** Phase variance and SNR before and after phase masking using PLS and PCM, evaluated at 5 Hz for phase-perturbed records under uniform and scaled amplitude conditions.
 
| Case | Method | Phase Variance (5 Hz) Before | Phase Variance (5 Hz) After | SNR Before | SNR After |
|:-----|:-------|:---:|:---:|:---:|:---:|
| Uniform amplitudes | PLS | 0.59 | 0.26 | -8.09 | -0.55 |
| | PCM | 0.59 | 0.26 | -8.09 | -0.56 |
| Scaled amplitudes | PLS | 0.59 | 0.32 | -9.14 | -3.60 |
| | PCM | 0.59 | 0.26 | -9.14 | -2.52 |
"""

# ╔═╡ 63b5c5fb-89aa-40ce-b6ed-7492a8a5e1b1
md"""
## 5. Conclusions
 
We introduced an amplitude-invariant phase-masking framework for coherence conditioning in scattered wavefields. The key contribution is phase masking via the circular mean (PCM), which estimates the dominant ensemble phase by normalizing spectra to unit magnitude prior to averaging, ensuring all records contribute equally regardless of their spectral energy.
 
Circular statistics provides a natural language for characterizing ensemble phase: the circular mean identifies the dominant phase direction across multichannel records, while circular variance quantifies the spread around it. Unlike stacking-based approaches, which estimate a representative phase but provide no measure of phase dispersion, this framework captures both.
 
Controlled synthetic experiments show that PCM and the conventional phase masking via local stack (PLS) are equivalent when amplitudes are uniform. When realistic trace-dependent amplitude variability is introduced, however, PLS becomes biased due to implicit amplitude weighting in the time-domain stack. PCM avoids this bias and consistently achieves lower phase variance across the bandwidth, producing more coherent summation.
 
These results demonstrate that amplitude-phase decoupling is helpful for reliable coherence conditioning in realistic multichannel data. Because the proposed formulation integrates directly with existing masking workflows, it provides a practical path toward improved phase stabilization and more robust coherence recovery in seismic imaging and other wave-based sensing applications.
"""

# ╔═╡ e65f1da4-d12b-4979-9b88-f27d9dc1706d
begin
    n_small = 100
    n_big   = 100_000
    freq_idx = 1:35
    f_axis   = positive_freq[freq_idx]

    function circ_stats_all_vs_firstN(phase::AbstractMatrix{<:Real}; n::Int=200)
        nf, nt = size(phase)
        nuse   = min(n, nt)
        μ_all  = zeros(Float64, nf)
        S_all  = zeros(Float64, nf)
        μ_n    = zeros(Float64, nf)
        S_n    = zeros(Float64, nf)
        for i in 1:nf
            g_all    = @view phase[i, :]
            g_n      = @view phase[i, 1:nuse]
            μ_all[i] = circ_mean(g_all)[1]
            S_all[i] = 1 - circ_r(g_all)
            μ_n[i]   = circ_mean(g_n)[1]
            S_n[i]   = 1 - circ_r(g_n)
        end
        return (; μ_all, S_all, μ_n, S_n, nuse)
    end

    # ── uniform amplitude stats 
    st_pert_u_100   = circ_stats_all_vs_firstN(phase_perturbed;        n=n_small)
    st_pert_u_big   = circ_stats_all_vs_firstN(phase_perturbed;        n=n_big)
    st_ls_u_100     = circ_stats_all_vs_firstN(predicted_phases_ls;    n=n_small)
    st_ls_u_big     = circ_stats_all_vs_firstN(predicted_phases_ls;    n=n_big)
    st_pcm_u_100    = circ_stats_all_vs_firstN(predicted_phases_pcm;   n=n_small)
    st_pcm_u_big    = circ_stats_all_vs_firstN(predicted_phases_pcm;   n=n_big)

    # ── variable amplitude stats 
    st_pert_s_100   = circ_stats_all_vs_firstN(phase_perturbed_scaled;        n=n_small)
    st_pert_s_big   = circ_stats_all_vs_firstN(phase_perturbed_scaled;        n=n_big)
    st_ls_s_100     = circ_stats_all_vs_firstN(predicted_phases_ls_scaled;    n=n_small)
    st_ls_s_big     = circ_stats_all_vs_firstN(predicted_phases_ls_scaled;    n=n_big)
    st_pcm_s_100    = circ_stats_all_vs_firstN(predicted_phases_pcm_scaled;   n=n_small)
    st_pcm_s_big    = circ_stats_all_vs_firstN(predicted_phases_pcm_scaled;   n=n_big)

    # ── 3a — uniform amplitudes 
    p3a = plot(f_axis, var_imposed[freq_idx],
               label="Input", lw=4, color=:black, linestyle=:dash,
               xlabel="Frequency (Hz)", ylabel="Phase Variance",
               title="(a) Phase perturbations with uniform amplitudes",
               grid=false, fontsize=12, tickfontsize=10,
               labelfontsize=12, legendfontsize=10)

    # perturbed
    plot!(p3a, f_axis, st_pert_u_big.S_all[freq_idx],
          label="V perturbed (100k)", lw=3, color=:red, linestyle=:solid)
    plot!(p3a, f_axis, st_pert_u_100.S_n[freq_idx],
          label="V perturbed (100)", lw=1, color=:red, linestyle=:dash)

    # PLS
    plot!(p3a, f_axis, st_ls_u_big.S_all[freq_idx],
          label="V local stack (100k)", lw=3, color=:blue, linestyle=:solid)
    plot!(p3a, f_axis, st_ls_u_100.S_n[freq_idx],
          label="V local stack (100)", lw=1, color=:blue, linestyle=:dash)

    # PCM
    plot!(p3a, f_axis, st_pcm_u_big.S_all[freq_idx],
          label="V circular mean (100k)", lw=3, color="#6fa76f", linestyle=:solid)
    plot!(p3a, f_axis, st_pcm_u_100.S_n[freq_idx],
          label="V circular mean (100)", lw=1, color="#6fa76f", linestyle=:dash)

    ylims!(p3a, 0.0, 0.8)

    # ── 3b — variable amplitudes 
    p3b = plot(f_axis, var_imposed[freq_idx],
               label="Input", lw=4, color=:black, linestyle=:dash,
               xlabel="Frequency (Hz)", ylabel="Phase Variance",
               title="(b) Phase perturbations with variable amplitudes",
               grid=false, fontsize=12, tickfontsize=10,
               labelfontsize=12, legendfontsize=10)

    # perturbed
    plot!(p3b, f_axis, st_pert_s_big.S_all[freq_idx],
          label="V perturbed (100k)", lw=3, color=:red, linestyle=:solid)
    plot!(p3b, f_axis, st_pert_s_100.S_n[freq_idx],
          label="V perturbed (100)", lw=1, color=:red, linestyle=:dash)

    # PLS
    plot!(p3b, f_axis, st_ls_s_big.S_all[freq_idx],
          label="V local stack (100k)", lw=3, color=:blue, linestyle=:solid)
    plot!(p3b, f_axis, st_ls_s_100.S_n[freq_idx],
          label="V local stack (100)", lw=1, color=:blue, linestyle=:dash)

    # PCM
    plot!(p3b, f_axis, st_pcm_s_big.S_all[freq_idx],
          label="V circular mean (100k)", lw=3, color="#6fa76f", linestyle=:solid)
    plot!(p3b, f_axis, st_pcm_s_100.S_n[freq_idx],
          label="V circular mean (100)", lw=1, color="#6fa76f", linestyle=:dash)

    ylims!(p3b, 0.0, 0.8)

    plot(p3a, p3b,
         layout=(1, 2),
         size=(1400, 500),
         left_margin=8Plots.mm,
         bottom_margin=8Plots.mm,
         top_margin=4Plots.mm,
         legend=false,
         dpi=300)
end

# ╔═╡ 345459f1-d417-4010-bdef-6a7dca1c6767
md"""
## Author Declarations
 
The authors have no conflicts of interest to disclose.
 
## Data Availability
 
The data that support the findings of this study were generated synthetically as described in the manuscript. The synthetic data and analysis scripts can be made available from the corresponding author upon reasonable request.
"""

# ╔═╡ 9d9bbc8a-ea28-465b-bd49-0ae338264d2b
md"""
## References
 
1. A. V. Oppenheim and J. S. Lim, "The importance of phase in signals," Proc. IEEE 69, 529–541 (1981).
2. H. L. Van Trees, Optimum Array Processing: Part IV of Detection, Estimation, and Modulation Theory (Wiley, New York, NY, 2002).
3. K. Paliwal, K. Wójcicki, and B. Shannon, "The importance of phase in speech enhancement," Speech Commun. 53, 465–494 (2011).
4. D. S. Williamson, Y. Wang, and D. Wang, "Complex ratio masking for monaural speech separation," IEEE/ACM Trans. Audio, Speech, Lang. Process. 24, 483–492 (2015).
5. R. M. Goldstein and C. L. Werner, "Radar interferogram filtering for geophysical applications," Geophys. Res. Lett. 25, 4035–4038 (1998).
6. A. Bakulin, I. Silvestrov, and D. Neklyudov, "Importance of phase guides from beamformed data for processing multi-channel data in highly scattering media," J. Acoust. Soc. Am. 147, EL447–EL452 (2020).
7. A. Bakulin, D. Neklyudov, and I. Silvestrov, "Seismic time-frequency masking for suppression of seismic speckle noise," Geophysics 88, V371–V385 (2023).
8. J. W. Goodman, Speckle Phenomena in Optics: Theory and Applications (Roberts and Company, Greenwood Village, CO, 2007).
9. A. Bakulin, D. Neklyudov, and I. Silvestrov, "Multiplicative random seismic noise caused by small-scale near-surface scattering and its transformation during stacking," Geophysics 87, V419–V435 (2022).
10. K. V. Mardia and P. E. Jupp, Directional Statistics (Wiley, New York, NY, 2009).
11. A. Rohatgi, A. Bakulin, and S. Fomel, "Data-driven analysis of seismic phase using circular statistics," The Leading Edge 44, 683–691 (2025).
12. A. Bakulin, I. Silvestrov, and M. Protasov, "Signal-to-noise ratio computation for challenging land single-sensor seismic data," Geophys. Prospect. 70, 629–638 (2022).
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CircStats = "2f6764a1-d620-4564-9394-76eb7c776766"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.10.15"
CircStats = "~1.0.4"
DataFrames = "~1.8.1"
Distributions = "~0.25.123"
FFTW = "~1.10.0"
HTTP = "~1.11.0"
Plots = "~1.41.5"
SpecialFunctions = "~2.7.2"
StatsBase = "~0.34.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "fa6dc16daadd6fe0ad553667673f0e965b26ee13"

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
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d0efe2c6fdcdaa1c161d206aa8b933788397ec71"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.6+0"

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
git-tree-sha1 = "c761b00e7755700f9cdf5b02039939d1359330e1"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.1.0"

[[deps.CommonSolve]]
git-tree-sha1 = "78ea4ddbcf9c241827e7035c3a03e2e456711470"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.6"

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
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

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
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

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
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "66381d7059b5f3f6162f28831854008040a4e905"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+1"

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
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

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
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

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
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

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
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

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
git-tree-sha1 = "b67985bd11331ccef26109a6269dbaae01474a72"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.11.6"

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
git-tree-sha1 = "67c6f1f085cb2671c93fe34244c9cccde30f7a26"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.5.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

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
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

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
git-tree-sha1 = "2ac022577e5eac7da040de17776d51bb770cd895"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

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
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

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
git-tree-sha1 = "1cc8ad0762e59e713ee3ef28f9b78b2c9f4ca078"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.5"

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
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

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
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

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
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "10a488dbecb88a9679c8f357d383d7d83dcc748d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.13"

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
git-tree-sha1 = "2700b235561b0335d5bef7097a111dc513b8655e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.2"

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
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"

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
git-tree-sha1 = "fa95b3b097bcef5845c142ea2e085f1b2591e92c"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.7.1"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
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
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

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
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

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
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

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
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

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

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

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
git-tree-sha1 = "45a20bd63e4fafc84ed9e4ac4ba15c8a7deff803"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.57+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

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
# ╟─27413a40-36df-11f1-3735-f9f362703981
# ╟─fd269c66-9602-4f90-81bb-ebb512c8bb35
# ╟─8f3256ed-c1b2-4e84-9763-af9ce49677d0
# ╟─e0caa217-e6f3-4013-bf62-b4f82e392c2c
# ╠═d6f47eee-dddf-4454-ac98-7fdf6f25ff99
# ╠═96f1b4dd-3e4c-493a-8b0e-57d11d7dea4a
# ╟─ad2b0861-50db-487f-b0fd-b6f95fa0358f
# ╟─3d3e9271-505e-42e9-9fac-52b54a02c8f4
# ╟─5073c1d5-9b7c-4636-8fa1-3c8882ac30b6
# ╟─19ea21f9-32b3-4e10-9439-d1731a0f5681
# ╟─5b2901db-9918-4bda-a084-0e4ada7e0301
# ╟─28465d38-af5e-4b70-b6b2-85785fc7c1cb
# ╟─74b5f360-988d-4f79-a96c-83682a02d2ae
# ╠═4153b364-dd01-4ddf-b556-4a2a592049d7
# ╠═40c12c7b-18be-4cf5-a42b-0ee57f30281b
# ╠═a2d7634f-9bcb-44cd-9f98-cb582cad0e3f
# ╠═fb919903-e6c4-45a4-97bd-b7e64b23c1ac
# ╠═a955bc18-5858-438f-b3d5-5821eceeb55a
# ╠═a0cb1efc-57b9-4fa2-aead-3633cb6405f3
# ╠═df2e462c-20dc-43a7-9b81-1be8245fcc96
# ╠═9324270e-ca58-46d9-a2ad-bbb7d259cf26
# ╠═3519b3d0-fa46-494b-9f38-70d6aeae0b82
# ╠═0173fbd2-044d-4eff-bbc0-e3c3ec67fa88
# ╠═44737414-985a-4cf5-80e5-ba822bb16231
# ╠═85e9ae48-a909-4c05-b8b5-d8d1ce39329f
# ╠═f0b9b033-5893-4e2e-8c72-d94547f50b11
# ╠═2b0bebc0-7159-4446-baec-ac5a5e5bb784
# ╠═3932ee58-1132-465a-9b03-33e690cecb81
# ╠═cebefee6-4653-48d9-893d-b221b0e83535
# ╠═31479529-c0cf-42f8-8fea-465084fe8c10
# ╠═07598432-1156-4d27-89c6-422dbaef6914
# ╠═ed98073b-8ac5-4bea-b0d6-cf646039196c
# ╠═5d44094f-b541-45f3-b164-a43af866d911
# ╠═dac97421-f42b-41e7-8569-a2cb5c5f2085
# ╠═f8b66581-f13a-451d-a77a-c19395d89ff1
# ╠═b3f6f24f-0260-4d0b-b80b-e60d5ffd2b3a
# ╟─3d6b96fd-53e9-4d89-8327-de4429cb83d8
# ╟─e91852af-5b48-41ca-b3ff-51d127d3d7e1
# ╟─c623aa1c-ebcf-4da4-a2c6-f417fb6abb75
# ╟─65e06a4c-287c-4a17-822b-a9d71de61368
# ╟─4b5f5011-156a-421e-ba81-e526f3d2a79e
# ╟─d4ae36d7-efa8-497e-b581-5268c0e3f863
# ╟─10eab8aa-247e-424d-81d3-e6feb0832d47
# ╟─8258d07d-4a86-4965-a509-e1645ca25ced
# ╟─f05c7e7f-4d99-444f-b179-d79375392896
# ╟─2cca3d97-5696-4fbe-ad34-368cbb4aa0b4
# ╠═65b1a5ba-c4e2-493d-a5cf-ec595fba83d7
# ╠═4a52c424-bc43-436c-9fc6-0bd4bb978a2a
# ╠═4ea2887e-b302-4377-bc4f-b838072576e4
# ╠═3e18f298-6ba0-472f-b421-c1bf043b5080
# ╟─46a4c8ea-33d0-46a8-aa63-4c3f1c04fad9
# ╟─4a115cf9-648f-4719-8d87-7d8e4796bdeb
# ╟─9bcaedf0-c58d-439e-933c-4761439bd1a4
# ╟─c755f57d-0b71-4e06-81ca-7b64e581d544
# ╠═36ea57c3-1ad5-473d-b2c1-15431754b670
# ╟─0573f916-bbb6-40e6-a6bd-8aea25ff36f1
# ╟─f8341612-b575-4e9c-bbdb-424c405664c8
# ╠═a4b60ad9-dd61-4980-bab7-373d7d6ba466
# ╠═b3bfcd66-2959-4192-a44b-4f27583ffe0d
# ╟─cdf2bd79-8064-488a-8eb5-d4690fc4e560
# ╠═07322a04-dfeb-4c12-97e8-b3e993baf352
# ╠═1948cc93-38ca-4be4-893a-531827b27312
# ╟─a506d20f-6b24-476c-9a46-8301c75c4b05
# ╟─63b5c5fb-89aa-40ce-b6ed-7492a8a5e1b1
# ╠═e65f1da4-d12b-4979-9b88-f27d9dc1706d
# ╟─345459f1-d417-4010-bdef-6a7dca1c6767
# ╟─9d9bbc8a-ea28-465b-bd49-0ae338264d2b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
