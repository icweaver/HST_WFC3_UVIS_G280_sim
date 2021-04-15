### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ e391b52a-92c5-4959-96f1-874dff1de951
begin
	using Revise
	import Pkg
	Pkg.activate("../Transits.jl/")
	Pkg.instantiate()
	using Plots, PlutoUI, Transits
end

# ╔═╡ 90217e93-4218-4463-8ac9-f5dab61bad1e
orbit_HP23 = KeplerianOrbit(
	0.99,              # ρₛ
	1.152,             # Rₛ
	2*1.2128867,         # P (fix bug)
	0.0,               # ecc
	0.0,               # t₀
	85.1 * π / 180.0,  # incl
)

# ╔═╡ b9c935f5-c842-407c-a8f8-b8f100b83b11
md"""
| Param      | Value | Definition
| ------------------ |:----------- |:-------
| ``N_\text{g}``      |  $(@bind N_g Slider(1:5, show_value=true)) | Number of groups
| ``N_\text{g}^{(\text{exp})}``      | $(@bind N_exp_g Slider(2:20, show_value=true)) | Number of exposures per group
| ``t_\text{start}``      |  $(@bind t_begin Slider(-3:0.5:-1, show_value=true)) hr | Time before mid-transit for first group
| ``t^{(\text{exp})}``   | $(@bind t_exp Slider(1:10:300, show_value=true)) s | Single exposure time
"""

# ╔═╡ 665df97b-bee1-4307-9c56-c0ac059fb575
md"""
An initial buffer of $(@bind t_buff Scrubbable(0:30, suffix=" minutes"))
"""

# ╔═╡ 7dd9c232-1784-427b-8bbb-627839b6cbb8
function compute_phases(t, P, t0)
	phase = ((t - t0) / P) % 1
	if phase ≥ 0.5
		phase = phase - 1.0
	end
	return phase
end

# ╔═╡ e82ad295-3505-44e8-822e-02730a27ac40
# https://discourse.julialang.org/t/findnearest-function/4143
function searchsortednearest(a, x)
   idx = searchsortedfirst(a,x)
   if (idx==1); return idx; end
   if (idx>length(a)); return length(a); end
   if (a[idx]==x); return idx; end
   if (abs(a[idx]-x) < abs(a[idx-1]-x))
      return idx
   else
      return idx-1
   end
end

# ╔═╡ 6404d888-9d97-11eb-318f-2f58cd9c2811
begin
	p_HP23 = 0.1113
	u_HP23 = [0.4, 0.26]
	ld_HP23 = PolynomialLimbDark(u_HP23)
	t_HP23 = range(0 - 3 / 24, 0 + 3 / 24, length=1000) # days
	ϕ = compute_phases.(t_HP23, orbit_HP23.P/2.0, orbit_HP23.t₀)
	fluxes_HP23 = @. ld_HP23(orbit_HP23, t_HP23, p_HP23)
	
	########
	# Groups
	########
	# In days
	P_HST = 0.0666667 # 96 min
	t_exp_s = t_exp / 3600 / 24 # 237 s
	
	########
	#N_g = 3
	#N_exp_g = 16
	t_group = N_exp_g * t_exp_s
	
	# Phase points sampled
	t = t_HP23
	t_start = t_begin / 24
	t_end = t_start + N_exp_g * t_exp_s
	t_group_1 = range(t_start, step=t_exp_s, length=N_exp_g)
	
	t_groups = [t_group_1]
	for i in 1:(N_g - 1)
		t_s = t_groups[i][end] + P_HST
	    t_e = t_s + t_group - t_exp_s
	    push!(t_groups, range(t_s, step=t_exp_s, length=N_exp_g))
	end
	
	phase_groups = []
	for t_group in t_groups
	    push!(phase_groups, compute_phases.(t_group, orbit_HP23.P/2.0, orbit_HP23.t₀))
	end
	
	phase_groups = hcat(phase_groups...)
	
	p = plot(xlabel="Planetary phase + 1", ylabel="Normalized flux")
	plot!(p, ϕ .+ 1, fluxes_HP23, label=false)
	scatter!(
		p,
		phase_groups .+ 1,
		#ones(size(phase_groups)[1]),
		fluxes_HP23[searchsortednearest.(Ref(t_HP23), hcat(t_groups...))],
		legend=false,
	)
end

# ╔═╡ a5517b00-6bce-468c-ab31-95476422f15e
APT_phase_range = compute_phases.([t_groups[begin][begin] - t_buff/60/24, t_groups[begin][begin+1] + t_buff/60/24], orbit_HP23.P/2.0, orbit_HP23.t₀) .+ 1;

# ╔═╡ 091c8ee7-ff12-4df4-97be-eb65bea3ae2e
md"""
would correspond to an APT phase of **$(APT_phase_range[1]) - $(APT_phase_range[2])** for the start of the first group 🚀
"""

# ╔═╡ Cell order:
# ╠═90217e93-4218-4463-8ac9-f5dab61bad1e
# ╟─b9c935f5-c842-407c-a8f8-b8f100b83b11
# ╟─6404d888-9d97-11eb-318f-2f58cd9c2811
# ╟─665df97b-bee1-4307-9c56-c0ac059fb575
# ╟─091c8ee7-ff12-4df4-97be-eb65bea3ae2e
# ╟─a5517b00-6bce-468c-ab31-95476422f15e
# ╟─7dd9c232-1784-427b-8bbb-627839b6cbb8
# ╟─e82ad295-3505-44e8-822e-02730a27ac40
# ╠═e391b52a-92c5-4959-96f1-874dff1de951
