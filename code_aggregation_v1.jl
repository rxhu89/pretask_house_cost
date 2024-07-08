"""
the aggregation scenarios
code created by Rongxing Hu

Julia v1.6.2
Atom v0.12.33
JuMP v0.21.8
Juno v0.8.4
"""
##
clearconsole() # only works in Juno, ctrl + D for clear variables
cd(dirname(@__FILE__)) # works in Atom or VSCode~ (but may not in Julia-self RE )

time_start = time()
##
import XLSX
using DataFrames
using Dates
using JuMP
using Cbc
using Plots
using Measures
using CSV
# using MAT
using Statistics

case_flag = "aggregated"

input_path = joinpath(pwd(), "Input_data" )
output_path = mkpath( joinpath(pwd(), "Output_data",case_flag) ) # if dir exists, return dir, otherwise, it creartes dir
## set parameters ================================================================
delta_t = 0.5; # ,time interval, in hour
N_hs = 3;     # the total number of houses
##-- component
soc_ini = 0.2 # the initial soc
soc_end = 0.2 # the end soc
soc_min = 0.1
soc_max = 1
eta_bat = 0.95 # the charging/discharging efficiency

Erate_bat1 = 4  # kWh
Prate_bat1 = 2  # kW

Erate_bat2 = 6 # kWh
Prate_bat2 = 3   # kW

Erate_bat3 = 8 # kWh
Prate_bat3 = 4 # kW

Erate_bat = [Erate_bat1, Erate_bat2, Erate_bat3]
Prate_bat = [Prate_bat1, Prate_bat2, Prate_bat3]

Prate_pv = [3, 0, 4] # no need to use since no Pgrid limits

## read input data ====================================================================
input_data_file = joinpath(input_path, "pre_task_data.xlsx")
input_data0 = XLSX.readxlsx(input_data_file)
input_data = input_data0[1][:]
# Pload_profile = DataFrame(input_data[3:end,2:4],:auto) # in kW
Pload_profile = input_data[3:end,2:4] # in kW
Ppv_profile = input_data[3:end,5:7]   # in kW
grid_price_imp = input_data[3:end,8:8] # in pence per kWh
grid_price_exp = input_data[3:end,9:9] # in pence per kWh

Pload_profile = Float64.(Pload_profile)
Ppv_profile  = Float64.(Ppv_profile)
grid_price_imp = Float64.(grid_price_imp)
grid_price_exp = Float64.(grid_price_exp)

T = size(grid_price_exp,1) # the total steps
M_big = 100*sum(maximum(eachrow(Pload_profile))); # set the big value

## modelling =====================================================================
m_opt = Model(Cbc.Optimizer)
# variables ==
@variable(m_opt, Pgrid_imp[1:T])   # imported power from grid, kW
@variable(m_opt, Pgrid_exp[1:T])   # exported power to grid, kW
@variable(m_opt, Pnet_house[1:T,1:N_hs]) # the net kw of each house
# battery
@variable(m_opt, E_bat[1:T,1:N_hs] )   # Battery energy
@variable(m_opt, Pch_bat[1:T,1:N_hs])  # charging power on grid side
@variable(m_opt, Pdis_bat[1:T,1:N_hs]) # discharging power on grid side
# PV
@variable(m_opt, Ppv[1:T,1:N_hs] )   # PV power

# objective expression, so can be adjusted later if we have more cases
objective_exp = @expression(m_opt, delta_t*sum( grid_price_imp[t]*Pgrid_imp[t] - grid_price_exp[t]*Pgrid_exp[t] for t = 1:T ))

# constraints- balance
@constraint(m_opt, pow_balan_limit[t=1:T,i=1:N_hs], Pnet_house[t,i] == Pload_profile[t,i]
                                                              - Pdis_bat[t,i] + Pch_bat[t,i] - Ppv[t,i]  )
@constraint(m_opt, pow_net_limit[t=1:T], sum(Pnet_house[t,i] for i = 1:N_hs) == Pgrid_imp[t] - Pgrid_exp[t] )
@constraint(m_opt, Pgrid_imp_limit[t=1:T], 0 <= Pgrid_imp[t] <= M_big )
@constraint(m_opt, Pgrid_exp_limit[t=1:T], 0 <= Pgrid_exp[t] <= M_big )
# constraints- battery
@constraint(m_opt, bat_Pch_limit[t=1:T,i=1:N_hs], 0 <= Pch_bat[t,i] <= Prate_bat[i] )
@constraint(m_opt, bat_Pdis_limit[t=1:T,i=1:N_hs], 0 <= Pdis_bat[t,i] <= Prate_bat[i] )

@constraint(m_opt, bat_energy_limit[t=1:T,i=1:N_hs], soc_min*Erate_bat[i] <= E_bat[t,i] <= soc_max*Erate_bat[i] )
@constraint(m_opt, bat_end_limit[t=T,i=1:N_hs], E_bat[t,i] == soc_ini*Erate_bat[i] );
@constraint(m_opt, bat_balan_limit_1[t=1,i=1:N_hs], E_bat[t,i] == soc_ini*Erate_bat[i] + (Pch_bat[t,i]*eta_bat - Pdis_bat[t,i]/eta_bat)*delta_t );
@constraint(m_opt, bat_balan_limit_2[t=2:T,i=1:N_hs], E_bat[t,i] == E_bat[t-1,i] + (Pch_bat[t,i]*eta_bat - Pdis_bat[t,i]/eta_bat)*delta_t );
# contraints -pv
@constraint(m_opt, pv_power_limit[t=1:T,i=1:N_hs], 0 <= Ppv[t,i] <= Ppv_profile[t,i] )

@objective(m_opt, Min, objective_exp);
## solve =================================================================
time_model = time()
time_model_execu = time_model - time_start
println("@@ Modelling finished in sec: $(time_model_execu)")
println(" **** Start solving the optimization ****** \n")

optimize!(m_opt);
sol_status = termination_status(m_opt);

time_sol_execu = time() - time_model
println("Optimization is solved in sec: $(time_sol_execu) \n")

if sol_status == MOI.INFEASIBLE
    println(" @@ Warning: ThE optimIzation solution is infeasible!");
end

if sol_status == MOI.OPTIMAL
    println( "(^^) Great: The optimization solution is feasible!")
    println("===== start Result Analysis ======= (^^)\n")

    sol_cost_total = JuMP.objective_value(m_opt)
    println("The objective function cost is \$", round(sol_cost_total; digits=3), ".\n");

    sol_Pgrid_imp = JuMP.value.(Pgrid_imp);
    sol_Pgrid_exp = JuMP.value.(Pgrid_exp);
    sol_Pnet_house = JuMP.value.(Pnet_house);
    sol_E_bat = JuMP.value.(E_bat);
    sol_Pch_bat = JuMP.value.(Pch_bat);
    sol_Pdis_bat = JuMP.value.(Pdis_bat);
    sol_Ppv = JuMP.value.(Ppv);

    sol_Pnet_pcc = sol_Pgrid_imp - sol_Pgrid_exp;
    sol_cost_step = sol_Pgrid_imp.*grid_price_imp*delta_t - sol_Pgrid_exp.*grid_price_exp*delta_t;

    sol_Pgrid_house_imp = zeros(T,N_hs);
    sol_Pgrid_house_exp = zeros(T,N_hs);
    for i = 1:N_hs
        for t = 1:T
            if sol_Pnet_house[t,i] >= 0
                sol_Pgrid_house_imp[t,i] = sol_Pnet_house[t,i]
            else
                sol_Pgrid_house_exp[t,i] = -1*sol_Pnet_house[t,i]
            end
        end
    end
    sol_cost_house_step = sol_Pgrid_house_imp.*grid_price_imp*delta_t - sol_Pgrid_house_exp.*grid_price_exp*delta_t

    # plot results ==========================================================
    vec_xaxis = collect(1:T)/2;
    # total
    plt1 = plot(vec_xaxis, sol_Pnet_pcc, label = "Grid Import", rightmargin=15mm)
    plot!(vec_xaxis, sum(Pload_profile, dims=2), label = "Total load" )
    plot!(vec_xaxis, sum(sol_Ppv, dims=2), label = "Total PV")
    plot!(vec_xaxis, sum(sol_Pch_bat-sol_Pdis_bat, dims=2), label = "Battery (ch+)", legend=:topleft)
    plot!(xlab="Hour", ylab="Power (kW)")
    plot!(twinx(),vec_xaxis,  -1*grid_price_exp, fillrange = grid_price_imp, fillalpha = 0.2, c = :gray,
                                     label = "Grid Price", ylab="Grid Price (Pence/kWh)", ylims=(-12,30) )
    display(plt1)
    savefig(joinpath(output_path,string("Total_overall","_",case_flag ,".png")) )
    # Battery E
    plt1 = plot(vec_xaxis, sum(sol_E_bat, dims=2), label = "Battery Energy", legend=:topleft, rightmargin=15mm)
    plot!(xlab="Hour", ylab="Battery Energy (kWh)")
    plot!(twinx(),vec_xaxis,  -1*grid_price_exp, fillrange = grid_price_imp, fillalpha = 0.2, c = :gray,
                                     label = "Grid Price", ylab="Grid Price (Pence/kWh)", ylims=(-12,30) )
    display(plt1)
    savefig(joinpath(output_path,string("Total_E","_",case_flag ,".png")) )

    # each grid net load and cost
    splt = plot(layout = grid(2, 1), rightmargin = 15mm )
    plot!(size=(600,600))
    plot!(splt[1], vec_xaxis, sol_Pnet_pcc, label = "Grid total" , rightmargin=15mm)
    plot!(splt[1], vec_xaxis, sol_Pnet_house, label = ["House 1" "House 2" "House 3"], legend=:topleft, ylims=(-3,9) )
    plot!(splt[1], title = "Net Load" )
    plot!(splt[1], xlab="Hour", ylab="Power (kW)")
    # plot!(ylim=(-3,9))
    plot!(twinx(splt[1]),vec_xaxis, -1*grid_price_exp, fillrange = grid_price_imp, fillalpha = 0.2, c = :gray,
                                     label = "Grid Price", ylabel = "Grid Price (Pence/kWh)", ylims=(-12,30) )

    plot!(splt[2], vec_xaxis, sol_cost_step, label = "Grid total", rightmargin=15mm)
    plot!(splt[2], vec_xaxis, sol_cost_house_step, label = ["House 1" "House 2" "House 3"], legend=:topleft, ylims=(-3*3,9*5) )
    plot!(splt[2], title = "Cost" )
    plot!(splt[2], xlab="Hour", ylab="Cost (Pence)")
     # plot!(ylim=(-3,9))
    plot!(twinx(splt[2]),vec_xaxis, -1*grid_price_exp, fillrange = grid_price_imp, fillalpha = 0.2, c = :gray,
                                                                      label = "Grid Price", ylabel = "Grid Price (Pence/kWh)", ylims=(-12,30) )
    plot!(background_color_legend=RGBA(1,1,1,0.6)) # 0.6 transparent
    display(splt)
    savefig(joinpath(output_path,string("each_Netload_cost","_",case_flag ,".png")) )

    # each house
    splt = plot(layout = grid(N_hs, 1), rightmargin = 15mm )
    plot!(size=(600,600))
    for i = 1:N_hs
        plot!(splt[i], vec_xaxis, sol_Pnet_house[:,i], label = string("Grid Import @ H ",i) )
        plot!(splt[i], vec_xaxis, Pload_profile[:,i], label = string("Load @ H ",i) )
        plot!(splt[i], vec_xaxis, sol_Ppv[:,i], label = string("PV @ H ",i))
        plot!(splt[i], vec_xaxis, sol_Pch_bat[:,i]-sol_Pdis_bat[:,i], label = string("Battery (ch+) @ H ",i), legend=:topleft)
        plot!(xlab="Hour", ylab="Power (kW)")
        plot!(splt[i], ylim = (-3, 4))
        plot!(twinx(splt[i]), vec_xaxis,  -1*grid_price_exp, fillrange = grid_price_imp, fillalpha = 0.2, c = :gray,
                                         label = "Grid Price", ylab="Grid Price (Pence/kWh)", ylims=(-12,30) )
    end
    plot!(background_color_legend=RGBA(1,1,1,0.6)) # 0.6 transparent
    display(splt)
    savefig(joinpath(output_path,string("each_overall","_",case_flag ,".png")) )


    # each battery
    splt = plot(layout = grid(N_hs, 1), rightmargin = 10mm )
    plot!(size=(600,600))
    for i = 1:N_hs
        plot!(splt[i], vec_xaxis,[sol_Pch_bat[:,i], -1*sol_Pdis_bat[:,i]], label = [string("charge @ H ", i) string("discharge @ H ", i)], legend=:topleft)
        plot!(splt[i], xlab = "Hour", ylab="Power (kW)")
        plot!(splt[i], ylim = (-maximum(Prate_bat), maximum(Prate_bat)))
        plot!(twinx(splt[i]), vec_xaxis, 0*sol_E_bat[:,i], fillrange = sol_E_bat[:,i], fillalpha = 0.2, c = :green, label = "E",
                                       ylabel = "Energy (kWh)", ylims=(0, maximum(Erate_bat)))
        # plot!(twinx(splt[i]), ylims = (0, maximum(Erate_bat)))
    end
    # plot!(title = "Battery Power" )
    plot!(background_color_legend=RGBA(1,1,1,0.6)) # 0.6 transparent
    display(splt)
    savefig(joinpath(output_path,string("each_Bat","_",case_flag ,".png")) )


    ## save data ====================================================================
    ## save csv from dataframe
    result_dispatch_df = DataFrame(hour = vec_xaxis, sol_Pnet_pcc = round.(sol_Pnet_pcc; digits=3),
                                              sol_Pgrid_imp = round.(sol_Pgrid_imp; digits=3),
                                              sol_Pgrid_exp = round.(sol_Pgrid_exp; digits=3) );
    for i = 1:N_hs
        result_pv_df = DataFrame(string("sol_Ppv", i) => sol_Ppv[:,i] );
        global result_dispatch_df = hcat(result_dispatch_df, round.(result_pv_df;digits=3))
    end
    for i = 1:N_hs
        result_bat_df = DataFrame(string("sol_E_bat", i) => sol_E_bat[:,i],
                                  string("sol_Pch_bat", i) => sol_Pch_bat[:,i],
                                  string("sol_Pdis_bat", i) => sol_Pdis_bat[:,i] );
        global result_dispatch_df = hcat(result_dispatch_df, round.(result_bat_df;digits=3))
    end
    for i = 1:N_hs
        result_Pnet_house_df = DataFrame(string("sol_Pnet_house", i) => sol_Pnet_house[:,i] );
        global result_dispatch_df = hcat(result_dispatch_df, round.(result_Pnet_house_df;digits=3))
    end
    for i = 1:N_hs
        sol_Pgrid_house_imp_df = DataFrame(string("sol_Pgrid_house_imp", i) => sol_Pgrid_house_imp[:,i] );
        global result_dispatch_df = hcat(result_dispatch_df, round.(sol_Pgrid_house_imp_df;digits=3))
    end
    for i = 1:N_hs
        sol_Pgrid_house_exp_df = DataFrame(string("sol_Pgrid_house_exp", i) => sol_Pgrid_house_exp[:,i] );
        global result_dispatch_df = hcat(result_dispatch_df, round.(sol_Pgrid_house_exp_df;digits=3))
    end
    for i = 1:N_hs
        sol_cost_house_step_df = DataFrame(string("sol_cost_house_step", i) => sol_cost_house_step[:,i] );
        global result_dispatch_df = hcat(result_dispatch_df, round.(sol_cost_house_step_df;digits=3))
    end
    vec_objective_Value = zeros(T);
    vec_objective_Value[1] =  sol_cost_total;
    vec_objective_Value_df = DataFrame(vec_objective_Value = round.(vec_objective_Value; digits=3));
    global result_dispatch_df = hcat(result_dispatch_df, round.(vec_objective_Value_df;digits=3))

    CSV.write( joinpath(output_path, string("result_dispatch","_",case_flag,".csv")), result_dispatch_df)


end
##
time_total = time() - time_start
println("***** Total Time Consumption in sec: $(time_total)")
