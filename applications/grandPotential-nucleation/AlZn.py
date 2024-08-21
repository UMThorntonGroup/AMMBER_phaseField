############################################################
import ammber.BinarySystems as BS
import matplotlib.pyplot as plt
from pycalphad import Database
plt.ion()
tdb_file = "Al_Li_Zn.tdb"
db = Database(tdb_file)
# pycalphad uses capitalized elements
elements = ["AL", "ZN"]
component = "ZN"
temperature = 570

AlZn_Sys = BS.BinaryIsothermalDiscreteSystem()
AlZn_Sys.fromTDB(db, elements, component, temperature)
AlZn_Fit = BS.BinaryIsothermal2ndOrderSystem()
AlZn_Fit.from_discrete_near_equilibrium(AlZn_Sys, x=0.35)

print(AlZn_Sys.phases.keys())
print(AlZn_Fit.phases.keys())
skew = AlZn_Sys.phases['FCC_A1'].xdata*0#8000
plt.plot(AlZn_Sys.phases['FCC_A1'].xdata,
         AlZn_Fit.phases['FCC_A1_0'].free_energy(AlZn_Sys.phases['FCC_A1'].xdata)+skew,
         linestyle='-', color='tab:orange',
         label="Fit 1",linewidth=2)

plt.plot(AlZn_Sys.phases['FCC_A1'].xdata,
         AlZn_Fit.phases['FCC_A1_1'].free_energy(AlZn_Sys.phases['FCC_A1'].xdata)+skew,
         linestyle='-', color='tab:blue',
         label="Fit 2",linewidth=2)

plt.plot(AlZn_Sys.phases['FCC_A1'].xdata,
         AlZn_Sys.phases['FCC_A1'].Gdata+skew,
         linestyle='-', color='k',
         label="FCC_A1",linewidth=1)

plt.xlim([0,1])
#plt.xticks([])
#plt.yticks([])

plt.savefig(fname="fit.png")

BS.write_binary_isothermal_parabolic_parameters(AlZn_Fit, output_file="AlZn.prm", template_file="template.prm", component="ZN")
