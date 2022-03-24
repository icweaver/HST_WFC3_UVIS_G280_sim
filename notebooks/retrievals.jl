### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 239a91a6-f68a-11eb-14fd-0ba8d08b08f9
begin
	import Pkg
	Pkg.activate(Base.current_project())

    using PlutoUI
	import MarkdownLiteral: @mdx
	using AlgebraOfGraphics, CairoMakie
	using CSV, DataFrames, DataFramesMeta, Glob, OrderedCollections
	using Statistics
	using PythonCall
end

# ╔═╡ 0132b4ab-0447-4546-b412-ec598b20d21d
@mdx """
# Retrievals
"""

# ╔═╡ a21aad0b-5998-420e-b437-d7ad262d0fe4
begin
	const DATA_DIR = "../data/retrievals"
	const FIG_DIR = "../figures/retrievals"
	TableOfContents()
end

# ╔═╡ 60dc161c-2aa2-4264-884d-6da3ead0e57b
@bind base_dir Select(glob("$(DATA_DIR)/*"))

# ╔═╡ d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
@bind fit_R0 Select(["NofitR0", "fitR0"])

# ╔═╡ 589afac8-0ea5-4962-b52b-7f035e91cf44
fname_suff = let
	b = basename(base_dir)
	label = occursin("offs", b) ? "offs" : "no_offs"
	label *= "_$(fit_R0)"
	label *= "_$(b)"
end

# ╔═╡ 093156c7-9da7-4814-9260-5173f27fa497
model_names = OrderedDict(
	"clear" => "NoHet_FitP0_NoClouds_NoHaze_$(fit_R0)",
	# "spot" => "Het_FitP0_NoClouds_NoHaze_$(fit_R0)",
	# "cloud" => "NoHet_FitP0_Clouds_NoHaze_$(fit_R0)",
	# "cloud+spot" => "Het_FitP0_Clouds_NoHaze_$(fit_R0)",
	# "haze" => "NoHet_FitP0_NoClouds_Haze_$(fit_R0)",
	#"clear+cloud+haze" => "NoHet_FitP0_Clouds_Haze_$(fit_R0)",
	# "haze+spot" => "Het_FitP0_NoClouds_Haze_$(fit_R0)",
	#"clear+spot+cloud+haze" => "Het_FitP0_Clouds_Haze_$(fit_R0)",
)

# ╔═╡ 0f65d095-09af-44d2-907b-c30e2c16b609
species = [
	# "Na",
	# "K",
	"TiO",
	# "Na_K",
	# "Na_TiO",
	# "K_TiO",
	# "Na_K_TiO",
	# "CO",
	# "H2O",
	# "NH3", # too low
	# "HCN", # too low
	# "CH4", # too low
	# "CO2", # too low
	# "FEH", # too low
]

# ╔═╡ a0094689-a9d5-4810-baba-bd7a96c27839
dashplus(x) = replace(x, '_' => '+')

# ╔═╡ 1c4fe72d-9872-4969-a62a-5163b5009bbb
@mdx """
## Retrieved transmission spectra 🐕
"""

# ╔═╡ 5569b57c-0585-4300-927b-5d089dde0f43
function plot_instr!(ax, fpath; color=:blue, label="add label")
	retr_instr = CSV.read(
		"$(base_dir)/$(fpath)", DataFrame;
		header = [:wav, :flux, :wav_d, :wav_u, :flux_err],
		comment = "#",
	)
	wav = retr_instr.wav .* 10_000
	errorbars!(ax, wav, retr_instr.flux, retr_instr.flux_err)
	scatter!(ax, wav, retr_instr.flux;
		markersize = 15,
		color,
		strokewidth = 1.5,
		label,
	)
end

# ╔═╡ db524678-9ee2-4934-b1bb-6a2f13bf0fa6
function get_retr_model(cube, sp, model)
	retr_model = cube[sp][model]["retr_model"]
	retr_model_sampled = cube[sp][model]["retr_model_sampled"]
	return retr_model, retr_model_sampled
end

# ╔═╡ cc011a66-37bd-4543-9a58-b11e1f785e52
function retrieval!(ax, model0, model_sampled; color=:blue, linewidth=3, label="")
	model = @rsubset model0 :wav < 1.0
	wav = model.wav .* 10_000
	wav_sampled = model_sampled.wav .* 10_000
	lines!(ax, wav, model.flux; color, linewidth, label)
	scatter!(ax, wav_sampled, model_sampled.flux;
		marker = :rect,
		markersize = 15,
		color,
		strokewidth = 1.5,
		strokecolor = color,
	)
	band!(ax, wav, model.flux_d, model.flux_u;
		color = (color, 0.25),
	)
end

# ╔═╡ 00a0f9c4-cd4d-4ae2-80b7-0c044239a571
function plot_retrieval!(ax, cube, sp, model; color=:blue, linewidth=3)
	retr_model, retr_model_sampled = get_retr_model(cube, sp, model)
	label = dashplus(sp) * " ($(model))"
	retrieval!(ax, retr_model, retr_model_sampled; color, linewidth, label)
end

# ╔═╡ 0f23e0d6-177d-4bf6-9660-f2c376b3146b
md"""
## Posterior distributions
"""

# ╔═╡ 1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
@mdx """
## Notebook setup 🔧
"""

# ╔═╡ eab74923-a084-468c-9b0d-c2cc98a23913
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig, pt_per_unit=1)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 44b3b8cd-4b83-4b27-a948-d1230489552f
begin
	@pyexec """
	global pickle, load_pickle
	import pickle

	def load_pickle(fpath):
		with open(fpath, "rb") as f:
			data = pickle.load(f)
		return data
	"""
	load_pickle(s) = @pyeval("load_pickle")(s)
end;

# ╔═╡ a7c68d25-a799-421b-9799-38837fa8a188
begin
	cube = OrderedDict()
	for sp ∈ species
		cube[sp] = OrderedDict()
		for (model_name, model_id) ∈ model_names
		cube[sp][model_name] = Dict()

		dirpath = "$(base_dir)/HATP23_E1_$(model_id)_$(sp)"

		cube[sp][model_name]["retr"] = load_pickle("$(dirpath)/retrieval.pkl")

		cube[sp][model_name]["retr_model"] = CSV.read(
			"$(dirpath)/retr_model.txt", DataFrame;
			header = [:wav, :flux, :flux_d, :flux_u],
			comment = "#",
		)

		cube[sp][model_name]["retr_model_sampled"] = CSV.read(
			"$(dirpath)/retr_model_sampled_Magellan_IMACS.txt", DataFrame;
			header = [:wav, :flux],
			comment = "#",
		)
		end
	end
	cube
end

# ╔═╡ 930ec094-7b11-48b8-818e-15c63ed6f8a5
dists = let
	data_py = cube["TiO"]["clear"]["retr"]["samples"]
	PyDict{String, Vector}(data_py)
end

# ╔═╡ 6ed6b481-4b84-485a-af2f-d2ecb943cd33
samples = pyconvert(Dict{String, Vector}, cube["TiO"]["clear"]["retr"]["samples"])

# ╔═╡ 0c477ad4-f3a7-4b72-af4f-5de693ddf2a9
begin
	dist = samples["logTiO"]
	hist(dist, bins=20)
	xlims!(-20, 0)
	current_figure()
end

# ╔═╡ e917bf4d-33b2-4ec2-b7ac-8a64c8a5f24a
mean(dist), std(dist)

# ╔═╡ e43f1834-73dd-4859-b847-f4c552561897
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (800, 600)
	const FIG_LARGE = (1_200, 1_000)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = let
		pal = Makie.ColorSchemes.Paired_8 |> reverse
		[pal[7:8] ; pal[5:6] ; pal[1:2]]
	end
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (
				xlabelsize = 18,
				ylabelsize = 18,
				topspinevisible = true,
				rightspinevisible = true,
				topspinecolor = :darkgrey,
				rightspinecolor = :darkgrey
			),
			Label = (
				textsize = 18,
				padding = (0, 10, 0, 0),
				font = AlgebraOfGraphics.firasans("Medium")
			),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			Text = (font = AlgebraOfGraphics.firasans("Medium"),),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)
	
	COLORS

end

# ╔═╡ e801501c-a882-4f2d-bbc1-40028c1c91d8
let
	fig = Figure(resolution=(1000, 500))
	ax = Axis(
		fig[1, 1],
		xlabel = "Wavelength (Å)",
		ylabel = "Transit depth (ppm)",
		#limits = (0.3, 1.1, 17_500, 21_000),
		limits = (2_000, 9_800, nothing, nothing),
		#limits = (4_600, 9_800, 17_000, 21_000)
		#limits = (4_000, 13_000, 15_500, 22_500),
		xticks = LinearTicks(7),
	)

	# for (i, sp) ∈ enumerate(species)
	# 	plot_retrieval!(ax, cube, sp, "clear"; color=COLORS[mod1(i, 6)], linewidth=1)
	# end
	# for (i, model) ∈ enumerate(models)
	# 	plot_retrieval!(ax, cube, "TiO", model; color=COLORS[i], linewidth=1)
	# end
	vlines!(ax, [5892.9, 7682.0, 8189.0], color=:grey, linestyle=:dash, linewidth=0.5)
	plot_retrieval!(ax, cube, "TiO", "clear"; color=COLORS[1], linewidth=3)
	
	fpath_suff = basename(base_dir)

	# Instrument
	plot_instr!(ax, "retr_JWST.txt";
		color = COLORS[5],
		label = "JWST",
	)
	plot_instr!(ax, "retr_Magellan_IMACS.txt";
		color = :white,
		label = "IMACS",
	)

	axislegend(orientation=:horizontal)

	savefig(fig, "$(FIG_DIR)/retr_$(fname_suff).pdf")
	
	fig
end

# ╔═╡ Cell order:
# ╟─0132b4ab-0447-4546-b412-ec598b20d21d
# ╠═a21aad0b-5998-420e-b437-d7ad262d0fe4
# ╟─60dc161c-2aa2-4264-884d-6da3ead0e57b
# ╟─d7ce97c1-82f2-46f1-a5ac-73e38e032fc8
# ╠═e917bf4d-33b2-4ec2-b7ac-8a64c8a5f24a
# ╠═0c477ad4-f3a7-4b72-af4f-5de693ddf2a9
# ╟─589afac8-0ea5-4962-b52b-7f035e91cf44
# ╟─093156c7-9da7-4814-9260-5173f27fa497
# ╠═0f65d095-09af-44d2-907b-c30e2c16b609
# ╠═a7c68d25-a799-421b-9799-38837fa8a188
# ╠═a0094689-a9d5-4810-baba-bd7a96c27839
# ╠═930ec094-7b11-48b8-818e-15c63ed6f8a5
# ╟─1c4fe72d-9872-4969-a62a-5163b5009bbb
# ╠═e801501c-a882-4f2d-bbc1-40028c1c91d8
# ╠═00a0f9c4-cd4d-4ae2-80b7-0c044239a571
# ╠═5569b57c-0585-4300-927b-5d089dde0f43
# ╠═db524678-9ee2-4934-b1bb-6a2f13bf0fa6
# ╠═cc011a66-37bd-4543-9a58-b11e1f785e52
# ╟─0f23e0d6-177d-4bf6-9660-f2c376b3146b
# ╠═6ed6b481-4b84-485a-af2f-d2ecb943cd33
# ╟─1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
# ╟─eab74923-a084-468c-9b0d-c2cc98a23913
# ╠═44b3b8cd-4b83-4b27-a948-d1230489552f
# ╠═e43f1834-73dd-4859-b847-f4c552561897
# ╠═239a91a6-f68a-11eb-14fd-0ba8d08b08f9
