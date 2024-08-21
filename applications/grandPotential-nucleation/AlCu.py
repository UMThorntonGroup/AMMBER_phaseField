import ammber.BinarySystems as bi
import matplotlib.pyplot as plt
from pycalphad import Database
tdb_file = "AlCu.TDB"
db = Database(tdb_file)
elements = ["AL", "CU"]
component = "CU"
temperature = 800

AlCu_Sys = bi.BinaryIsothermalDiscreteSystem()
AlCu_Sys.fromTDB(db, elements, component, temperature)
AlCu_Sys_neareq = AlCu_Sys.resample_near_equilibrium(0.1)


AlCu_Sys_neareq.phases['LIQUID'] = AlCu_Sys.phases['LIQUID'].resample_near_xpoint(0.1)
AlCu_Fit = bi.BinaryIsothermal2ndOrderSystem()
AlCu_Fit.from_discrete(AlCu_Sys_neareq)

print(AlCu_Fit.phases.keys())
x = AlCu_Sys.phases['FCC_A1'].xdata
plt.plot(x, AlCu_Fit.phases['FCC_A1_0'].free_energy(x),
            linestyle=':', color='tab:orange', label="FCC_A1_0 fit",linewidth=1)

plt.plot(AlCu_Sys.phases['FCC_A1'].xdata, AlCu_Sys.phases['FCC_A1'].Gdata,
            linestyle='-', color='tab:orange',label="FCC_A1",linewidth=1)

plt.plot(x, AlCu_Fit.phases['AL2CU_1'].free_energy(x),
            linestyle=':', color='tab:green',linewidth=1, label="Al2Cu_C16 fit")

plt.plot(AlCu_Sys.phases['AL2CU'].xdata, AlCu_Sys.phases['AL2CU'].Gdata,
            linestyle='-', color='tab:green',linewidth=1, label="Al2Cu")

plt.plot(x, AlCu_Fit.phases['LIQUID'].free_energy(x),
            linestyle=':', color='tab:blue',linewidth=1, label="Liquid fit")

plt.plot(AlCu_Sys.phases['LIQUID'].xdata, AlCu_Sys.phases['LIQUID'].Gdata,
            linestyle='-', color='tab:blue',linewidth=1, label="Liquid")

plt.xlim([0.0, 0.4])
plt.ylim([-50000.0, -25000])
plt.xlabel("X_Cu")
plt.ylabel("G (J/mol)")
plt.savefig("AlCu_eutectic.png")

bi.write_binary_isothermal_parabolic_parameters(AlCu_Fit, output_file="AlCu.prm", 
                                                template_file="template.prm", phases=["LIQUID", "FCC_A1_0", "AL2CU_1"], component="CU")