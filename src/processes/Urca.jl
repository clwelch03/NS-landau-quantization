module Urca
export rate_wrapper, emissivity_plot, rb_plot
include("DirectUrca.jl")
include("ModifiedUrca.jl")
include("../utils/Constants.jl")
include("../eos/EOS.jl")
include("../utils/HelperFunctions.jl")
using .DirectUrca, .ModifiedUrca, .EOS, .HelperFunctions, PyPlot, Printf, Optim

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 12

LINESTYLES = ["dotted", "dashed", "dashdot"]
ST_GREEN = "#3DAE2B"
ST_BLUE = "#00A0DF"
ST_PINK = "#ED40A9"
ST_ORANGE = "#F38B00"
ST_PURPLE = "#B14FC5"

COLORS = [ST_GREEN, ST_BLUE, ST_PINK, ST_ORANGE, ST_PURPLE]


"""
    rate_wrapper(mu_B, B, T; mageos=true, nuc_inter=true, spin_split=true, thermal_population=true, analytical_approx=false, Mred_in_int=true)

Compute the direct Urca, modified Urca, and quasiclassical DUrca rates, as well as the density and number of available LLs, for a given baryon chemical potential.
"""
function rate_wrapper(mu_B, B, T; mageos=true, nuc_inter=true, spin_split=true, thermal_population=true, analytical_approx=false, Mred_in_int=true)
    # calculate EoS params
    if mageos
        eos_data = eos_data_mag(iufsu_star_constants_mag, mu_B/197.3, B/1e15) * 197.3
    else
        eos_data = eos_data_rmf(iufsu_star_constants, mu_B/197.3) * 197.3
    end

    # process EoS params 
    if nuc_inter
        k_Fn, k_Fp, µ_e, M_n_star, M_p_star = process_eos_data(eos_data, mageos)
        E_Fn = mu_B
        E_Fp = mu_B - µ_e 
    else
        k_Fn, k_Fp, µ_e, _, _ = process_eos_data(eos_data, mageos)
        M_n_star = 940.6
        M_p_star = 938.3
        E_Fn = mu_B
        E_Fp = mu_B - µ_e
    end
    density_val = density(k_Fn, k_Fp)

    # calculate rate values for both reactions
    DUrca_rate_val, n_maxes = direct_urca_rate(mu_B, k_Fn, k_Fp, µ_e, M_n_star, M_p_star, B, T;
        spin_split=spin_split, thermal_population=thermal_population, analytical_approx=analytical_approx, Mred_in_int=Mred_in_int, cgs=true)
    MUrca_rate_val = modified_urca_rate(k_Fn, k_Fp, M_n_star, M_p_star, T, cgs=true)
    qc_rate_val = qc_emissivity(k_Fn, k_Fp, M_n_star, M_p_star, B, T, µ_e=µ_e)[1]

    # println("durca rate args: ", mu_B, " ", k_Fn, " ", k_Fp, " ", µ_e, " ", M_n_star, " ", M_p_star)
    DUrca_rate_val, MUrca_rate_val, qc_rate_val, density_val, n_maxes #, durca_threshold_val
end


"""
    emissivity_plot(B, T; mu_range=[950, 980], points=100,
                                            include_MUrca=true,
                                            include_BY=false,
                                            include_new_LL=false,
                                            log_plot=true,
                                            mageos=true,
                                            nuc_inter=true,
                                            spin_split=true,
                                            thermal_population=true, 
                                            analytical_approx=false,
                                            Mred_in_int=true,
                                            refine_points=true,
                                            debug_title_info=false,
                                            return_data=true)

Create a plot comparing emissivity at variable density and fixed B.
"""
function emissivity_plot(B, T; mu_range=[950, 980], points=100,
                        include_MUrca=true,
                        include_BY=false,
                        include_new_LL=false,
                        log_plot=true,
                        mageos=true,
                        nuc_inter=true,
                        spin_split=true,
                        thermal_population=true, 
                        analytical_approx=false,
                        Mred_in_int=true,
                        refine_points=true,
                        debug_title_info=false,
                        return_data=true)
    
    B_str = @sprintf("%.2e", B) # gives us scientific notation

    if debug_title_info
        title_text = "Emissivity vs density: B = $B_str G, T = $T MeV\nSpin split $spin_split, therm pop $thermal_population, mageos $mageos, nuc. inter. $nuc_inter"
    else
        title_text = "Emissivity vs density: B = $B_str G, T = $T MeV"
    end

    # initialize arrays to store rate values and other info
    mu_vals = range(0, stop=(mu_range[2]-mu_range[1])^(5/9), length=points).^1.8 .+ mu_range[1]
    DUrca_rate_vals = Array{Float64}(undef, points)
    MUrca_rate_vals = Array{Float64}(undef, points)
    qc_rate_vals = Array{Float64}(undef, points)
    density_vals = Array{Float64}(undef, points)
    if spin_split # n_p_up, n_p_down, n_e
        old_n_maxes_arr = Array{Int16}(undef, (points, 3))
        n_maxes_arr = Array{Int16}(undef, (points, 3))
        old_n_maxes_arr[1, :] = [-1 -1 -1]
    else # n_p, n_e
        old_n_maxes_arr = Array{Int16}(undef, (points, 2))
        n_maxes_arr = Array{Int16}(undef, (points, 2))
        old_n_maxes_arr[1, :] = [-1 -1]
    end

    durca_wrapper(mu_B) = -rate_wrapper(mu_B, B, T; 
                                            mageos=mageos, nuc_inter=nuc_inter, spin_split=spin_split,
                                            analytical_approx=analytical_approx, Mred_in_int=false)[1] # for optimization

    # compute rate values at initial mu_B points
    @debug "Computing initial grid of rate values"
    for point in 1:points
        if point % 10 == 0
            println("Computing point $point of $points \r")
        end
        durca_val, murca_val, qc_val, density_val, n_max_val = rate_wrapper(mu_vals[point], B, T; mageos=mageos, nuc_inter=nuc_inter, spin_split=spin_split,
                                thermal_population=thermal_population, analytical_approx=analytical_approx, Mred_in_int=Mred_in_int)

        # update DUrca, MUrca, BY, and density arrays
        DUrca_rate_vals[point] = durca_val
        MUrca_rate_vals[point] = murca_val
        qc_rate_vals[point] = qc_val
        density_vals[point] = density_val

        # put current n_max in arrays as n_max and 'next' old_n_max
        n_maxes_arr[point, :] = n_max_val
        if point != points
            old_n_maxes_arr[point+1, :] = n_max_val
        end
    end

    # find points *close* to maxima and put them *at* the maxima
    if refine_points
        println("Refining points near maxima")
        for point in 2:points-1
            if DUrca_rate_vals[point] > DUrca_rate_vals[point-1] && DUrca_rate_vals[point] > DUrca_rate_vals[point+1]
                new_mu = Optim.minimizer(optimize(durca_wrapper, mu_vals[point-1], mu_vals[point+1]))
                new_durca_val, new_murca_val, new_qc_val, new_density_val, new_n_max_val = rate_wrapper(new_mu, B, T; mageos=mageos, nuc_inter=nuc_inter, spin_split=spin_split,
                                    thermal_population=thermal_population, analytical_approx=analytical_approx, Mred_in_int=Mred_in_int)

                # update DUrca, MUrca, BY, and density arrays
                DUrca_rate_vals[point] = new_durca_val
                MUrca_rate_vals[point] = new_murca_val
                qc_rate_vals[point] = new_qc_val
                density_vals[point] = new_density_val

                # put current n_max in arrays as n_max and 'next' old_n_max
                n_maxes_arr[point, :] = new_n_max_val
                old_n_maxes_arr[point+1, :] = new_n_max_val
            end
        end
    end

    # plot emissivity curves
    plot(density_vals, DUrca_rate_vals, linewidth=2.5, linestyle="-", color=ST_GREEN, label="DUrca rate")
    if include_MUrca
        plot(density_vals, MUrca_rate_vals, linewidth=2.5, linestyle="--", color=ST_BLUE, label="MUrca rate")
    end
    if include_BY
        plot(density_vals, qc_rate_vals, linewidth=2.5, linestyle=":", color=ST_PINK, label="BY DUrca rate")
    end

    # show where new LLs become available
    if include_new_LL
        new_n_mask = old_n_maxes_arr .!= n_maxes_arr
        new_n_mask[1, :] = spin_split ? [0 0 0] : [0 0] # first point technically has new n values but like. cmon
        y_min, y_max = ylim()
        if spin_split
            new_n_p_up_mask = new_n_mask[:, 1]
            new_n_p_down_mask = new_n_mask[:, 2]
            new_n_e_mask = new_n_mask[:, 3]
            vlines(density_vals[new_n_p_up_mask], ymin=y_min, ymax=y_max, linewidth=1, linestyle=":", color="gray", label="New LL, p spin up")
            vlines(density_vals[new_n_p_down_mask], ymin=y_min, ymax=y_max, linewidth=1, linestyle="-.", color="gray", label="New LL, p spin down")
            vlines(density_vals[new_n_e_mask], ymin=y_min, ymax=y_max, linewidth=1, linestyle="--", color="gray", label="New LL, e")
        else
            new_n_p_mask = new_n_mask[:, 1]
            new_n_e_mask = new_n_mask[:, 2]
            vlines(density_vals[new_n_p_mask], ymin=y_min, ymax=y_max, linewidth=1, linestyle=":", color="gray", label="New LL, p")
            vlines(density_vals[new_n_e_mask], ymin=y_min, ymax=y_max, linewidth=1, linestyle="--", color="gray", label="New LL, e")
        end
    end

    # make the plot look good
    title(title_text)
    if log_plot
        plt.gca()[:set_yscale]("log")
    end
    xlabel(L"Density (fm$^{-3}$)")
    ylabel(L"Emissivity (erg cm$^{-3}$ s$^{-1}$)")
    legend()
    show()
    PyPlot.display_figs()

    if return_data
        density_vals, DUrca_rate_vals, MUrca_rate_vals, qc_rate_vals, n_maxes_arr
    end
end


"""
    rb_plot(x_range, B, temp, points; toys=false, analytical_approx=false, Mred_in_int=false)

Create a plot comparing our R_B (at various `temps`) to the R_B from Baiko and Yakovlev (1999).

Can be run with or without our extra features (e.g. spin splitting) according to the `toys` toggle.
"""
function rb_plot(x_range, B, temp, points; toys=false, analytical_approx=false, Mred_in_int=false)
    x_vals = Array{Float64}(undef, points)
    µ_bounds = find_µ_range_endpoints(x_range[1], x_range[2], B)
    µ_vals = range(0, stop=(µ_bounds[2]-µ_bounds[1])^(5/9), length=points).^1.8 .+ µ_bounds[1]

    qc_RB_vals = Array{Float64}(undef, points)
    our_RB_vals = Array{Float64}(undef, points)
    density_vals = Array{Float64}(undef, points)

    # compute values at initial points
    @debug "Computing initial grid of points"
    Threads.@threads for point in 1:points
        our_RB_vals[point], qc_RB_vals[point], x_vals[point], density_vals[point] = rb(µ_vals[point], B, temp; toys=toys, analytical_approx=analytical_approx, Mred_in_int=Mred_in_int)
        print(point, ", ")
    end

    print("\n")
    println(x_vals)
    println(density_vals)
    println(qc_RB_vals)
    println(our_RB_vals)
end

end # module Urca