## computes the values of eta per region and sector based on first iteration
## values of the labor demand and parameters from excel sheet

print("Using labor allocation enn's page 36 2021 feb version:")

labor = {}
labor["r1"] = {}
labor["r1"][0] = 0.302619997111055
labor["r1"][1] = 0.00419116881427591
labor["r1"][2] = 0.00252919056260904
labor["r1"][3] = 0.000659643512059885
labor["r2"] = {}
labor["r2"][0] = 0.38079407548916
labor["r2"][1] = 0.00324868632569556	
labor["r2"][2] = 0.00170864645984442	
labor["r2"][3] = 0.0040785917252996
labor["r3"] = {}
labor["r3"][0] = 0.294858949486565
labor["r3"][1] = 0.0040377664987181
labor["r3"][2] = 0.0016112222065384
labor["r3"][3] = 0.00108206180817834
labor["r4"] = {}
labor["r4"][0] = 0.322415582463319
labor["r4"][1] = 0.00198363305125737
labor["r4"][2] = 0.010418884020606
labor["r4"][3] = 0.000431900464817843
labor["r5"] = {}
labor["r5"][0] = 0.380002065054839
labor["r5"][1] = 0.00565573224760838
labor["r5"][2] = 0.00223812439446772
labor["r5"][3] = 0.000744078303084591
labor["r6"] = {}
labor["r6"][0] = 0.400724363509264
labor["r6"][1] = 0.00313252655162304
labor["r6"][2] = 0.00317351309584523
labor["r6"][3] = 0.00404459684326768

alpha = [0.3               , 0.272           , 0.391            , 0.82]
nu    = [0.0812233567326101, 0.58730958617077, 0.205526512322629, 0   ]

eta1 = {}
for r in ["r1","r2","r3","r4","r5","r6"]:
	eta1[r] = {}
	for x in range(1,4):
		eta1[r][x] = (1-alpha[0]-nu[0])/(1-alpha[x]-nu[x]) * (labor[r][x]/(labor[r][0]*nu[0]))
		print("eta",r,x,eta1[r][x])


print("-"*15)
print("Using resource consumption (LF version of eqn for X, page 37) in 2021 feb version:")

ems = {}
ems["r1"] = {}
ems["r1"][1] = 13.9
ems["r1"][2] = 9.6
ems["r2"] = {}
ems["r2"][1] = 10.3
ems["r2"][2] = 6.2
ems["r3"] = {}
ems["r3"][1] = 12.7
ems["r3"][2] = 5.8
ems["r4"] = {}
ems["r4"][1] = 5.3
ems["r4"][2] = 31.86
ems["r5"] = {}
ems["r5"][1] = 15.5
ems["r5"][2] = 7.02
ems["r6"] = {}
ems["r6"][1] = 7.4
ems["r6"][2] = 8.58

ys = {}
ys["r1"] = 0.1567
ys["r2"] = 0.1885
ys["r3"] = 0.1448
ys["r4"] = 0.1345
ys["r5"] = 0.1626
ys["r6"] = 0.1478

evv = [0.0004033, 0.000043]

eta2 = {}
for r in ["r1","r2","r3","r4","r5","r6"]:
	s = 0
	eta2[r] = {}
	for x in range(1,3):
		eta2[r][x] = evv[x-1]*ems[r][x] / (nu[0]*nu[x]*ys[r])
		s         += eta2[r][x]
		print("eta",r,x,eta2[r][x])
	eta2[r][3] = 1-s
	print("eta",r,3,eta2[r][3])


print("-"*15)
print("Using energy ratio (equation A.3) from 2021 feb version:")

kap = [0.533249934081155, 0.29482640371216, 0.171923662207]
rho = 0.5 

ems = {}
ems["r1"] = {}
ems["r1"][1] = 0.550677508356912
ems["r1"][2] = 0.0797660371305951
ems["r1"][3] = 0.0801714870104502
ems["r1"][0] = sum([kap[x-1]*ems["r1"][x]**rho for x in range(1,4)])**(1./rho)
ems["r2"] = {}
ems["r2"][1] = 0.166806430900318
ems["r2"][2] = 0.0183539494674005
ems["r2"][3] = 1.54522602611472
ems["r2"][0] = sum([kap[x-1]*ems["r2"][x]**rho for x in range(1,4)])**(1./rho)
ems["r3"] = {}
ems["r3"][1] = 0.310653906972052
ems["r3"][2] = 0.0196758352386005
ems["r3"][3] = 0.131121187605852
ems["r3"][0] = sum([kap[x-1]*ems["r3"][x]**rho for x in range(1,4)])**(1./rho)
ems["r4"] = {}
ems["r4"][1] = 0.0179755416969855
ems["r4"][2] = 0.197255718599036
ems["r4"][3] = 0.00500842270236741
ems["r4"][0] = sum([kap[x-1]*ems["r4"][x]**rho for x in range(1,4)])**(1./rho)
ems["r5"] = {}
ems["r5"][1] = 0.165189871960399
ems["r5"][2] = 0.0102896844480273
ems["r5"][3] = 0.0168041803505749
ems["r5"][0] = sum([kap[x-1]*ems["r5"][x]**rho for x in range(1,4)])**(1./rho)
ems["r6"] = {}
ems["r6"][1] = 0.0143061508847038
ems["r6"][2] = 0.00584038866566252	
ems["r6"][3] = 0.140170968547232
ems["r6"][0] = sum([kap[x-1]*ems["r6"][x]**rho for x in range(1,4)])**(1./rho)

eta3 = {}
for r in ["r1","r2","r3","r4","r5","r6"]:
	eta3[r] = {}
	for x in range(1,4):
		eta3[r][x] = kap[x-1] * (ems[r][x]/ems[r][0])**rho
		print("eta",r,x,eta3[r][x])

print("-"*15)
print("Compare three etas for region=1 sector=1:")
print("labor alloc =>",eta1["r1"][1])
print("res consump =>",eta2["r1"][1])
print("enrgy ratio =>",eta3["r1"][1])


