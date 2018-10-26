
import string
import numpy as np
import sympy
import csv

def SCC(s): #strip control codes
	return ''.join([x for x in s if x in string.printable])

def readCSV(fileName,delimiter=","):
	dataDict = {}
	with open(fileName+".csv") as air_Cp_file:
		csvReader = csv.reader( air_Cp_file, delimiter="," )
		currentLine = 0
		keys = []
		for row in csvReader:
			if (currentLine==0):
				keys = [ SCC(s) for s in row]
				for key in keys:
					dataDict[key] = []
			else:
				for i in range(len(row)):
					dataDict[keys[i]].append(float(row[i]))

			currentLine+=1

	return dataDict

def between(x,a,b):
	if (x>=min(a,b) and x<=max(a,b)):
		return True
	return False

def interpolate(x,x1,x2,y1,y2): #interpolate linear
	p = (x-x1)/(x2-x1)
	return p*(y2-y1)+y1

def interpolateFromData(x,xData,yData): #interpolate a point on a data
	#the input data must be ordered
	for i in range(len(xData))[1:]:
		if between(x,xData[i-1],xData[i]):
			return interpolate(x,xData[i-1],xData[i],yData[i-1],yData[i])
	print("x is out of bounds of the function!")
	return None

air_data = readCSV("air_properties")

def get_air_Cp(T): #T in K
	T_data = air_data["T"]
	Cp_data = air_data["Cp"]
	return interpolateFromData(T,T_data,Cp_data)

def get_air_Pr(T):
	T_data = air_data["T"]
	Pr_data = air_data["Pr"]
	return interpolateFromData(T,T_data,Pr_data)

def get_air_Rho(T): 
	T_data = air_data["T"]
	Rho_data = air_data["rho"]
	return interpolateFromData(T,T_data,Rho_data)

def get_air_u(T): #dynamic viscocity
	T_data = air_data["T"]
	u_data = air_data["u"]
	return interpolateFromData(T,T_data,u_data)

jH_andCf_data = readCSV("jH_and_Cf")

def get_jH_and_Cf(Re):
	Re_data = jH_andCf_data["Re"]
	jH_data = jH_andCf_data["jH"]
	Cf_data = jH_andCf_data["Cf"]
	jH = interpolateFromData(Re,Re_data,jH_data)
	Cf = interpolateFromData(Re,Re_data,Cf_data)
	return (jH,Cf)

def C_to_K(C):
	return C+273.15

def calc_C(massFlow,specificHeat):
	return massFlow*specificHeat

def calc_Cr(C_h,C_c):
	C_min = float(min(C_h,C_c))
	C_max = float(max(C_h,C_c))
	return C_min/C_max

def logMeanDeltaT(T_max_m,T_min_i,T_min_o):
	pass

def e_crossflow_both_unmixed(NTU,Cr):
	return 1.0 - np.exp( (1.0/Cr)*np.power(NTU,0.22)*( np.exp(-Cr*np.power(NTU,0.78)) - 1.0 ) )

def calculate_G(rho_m,Pr,dP,ntu,jH,Cf,n=0.8):
	return np.sqrt( ( 2.0/( (1.0/rho_m)*(Pr**(2.0/3.0)) ) ) * ( (n*dP)/ntu ) * (jH/Cf) )

def approx_rho_m(rho_i,rho_o):
	return 1.0 / ( 0.5*( 1.0/rho_i + 1.0/rho_o ) )

epsilon = 0.8381

T_g_i = C_to_K(900.0) #k
m_dot_g = 1.66 #kg/s
P_g_i = 160e3 #Pa
dP_g = 9.05e3 #Pa
P_g_o = P_g_i - dP_g

T_a_i = C_to_K(200) #k
m_dot_a = 2.0 #kg/s
P_a_i = 200e3 #Pa
dP_a = 8.79e3 #Pa
P_a_o = P_a_i - dP_a

#wall parameters (the plate between the flows)
wall_k = 18.0 #thermal conductivity
wallThickness = 0.5e-3 #m

#fin parameters
fin_k = 18.0 #thermal conductivity of the material
D_h = 1.54e-3 #[m] hydraulic diameter
r_h = D_h/4.0
plateSpacing = 2.49e-3 #m
finThickness = 0.102e-3 #m
finLength = 3.18e-3 #m
finHeight = plateSpacing/2.0 - finThickness
wettedPerimeter = 2.0*(finLength+finThickness)
A_k = finLength*finThickness
finAreaRatio = 0.785
alpha = 941.6 #m^2/m^3

#exit pressure loss coefficients
K_c = 0.36
K_e = 0.42

C_r = min(m_dot_g,m_dot_a)/max(m_dot_g,m_dot_a) #initial approximation (assuming Cp's are equivalent)

T_g_o = T_g_i - epsilon*(T_g_i-T_a_i)
T_a_o = T_a_i + epsilon*(T_g_i-T_a_i)*C_r

T_g_m = (T_g_i+T_g_o)/2.0 #linear (very crude) approximation
T_a_m = (T_a_i+T_a_o)/2.0

C_p_g = 1e3*get_air_Cp(T_g_m) #convert from Kj/kg
C_p_a = 1e3*get_air_Cp(T_a_m)

Pr_g = get_air_Pr(T_g_m)
Pr_a = get_air_Pr(T_a_m)

Rho_g = get_air_Rho(T_g_m)
Rho_a = get_air_Rho(T_a_m)
Rho_g_i = get_air_Rho(T_g_i)
Rho_a_i = get_air_Rho(T_a_i)
Rho_g_o = get_air_Rho(T_g_o)
Rho_a_o = get_air_Rho(T_a_o)

u_g = 1e-5*get_air_u(T_g_m) #table data in: 10^-5 kg/ms, thus need to convert
u_a = 1e-5*get_air_u(T_a_m)

C_g = m_dot_g*C_p_g
C_a = m_dot_a*C_p_a

C_min = min(C_g,C_a)
C_max = max(C_g,C_a)

C_r = C_min/C_max

q = epsilon*C_min*(T_g_i-T_a_i)

NTU_guess = -np.log(1-epsilon) #initial guess (based on heat exchanger with Cr=0)
NTU_sym = sympy.Symbol("NTU_sym")
f = 1.0 - sympy.exp( (1.0/C_r)*(NTU_sym**0.22)*( sympy.exp(-C_r*(NTU_sym**0.78)) - 1.0 ) ) - epsilon
NTU = float(sympy.nsolve(f, NTU_guess))

#for gas heat-exchangers: ntu_h = ntu_c = 2*NTU
ntu_g = ntu_a = 2.0*NTU

jH_over_Cf_guess = 0.3
jH_g = jH_over_Cf_guess
Cf_g = 1.0
jH_a = jH_over_Cf_guess
Cf_a = 1.0

G_g = None
G_a = None
for i in range(100):
	Rho_g_m = approx_rho_m(Rho_g_i,Rho_g_o)
	Rho_a_m = approx_rho_m(Rho_a_i,Rho_a_o)
	G_g = calculate_G(Rho_g_m,Pr_g,dP_g,ntu_g,jH_g,Cf_g)
	G_a = calculate_G(Rho_a_m,Pr_a,dP_a,ntu_a,jH_a,Cf_a)
	Re_g = (G_g*D_h)/u_g
	Re_a = (G_a*D_h)/u_a
	jH_g,Cf_g = get_jH_and_Cf(Re_g)
	jH_a,Cf_a = get_jH_and_Cf(Re_a)

L_g = None
L_a = None
L = None

iterations = 10
for i in range(iterations):
	print("========Interation:"+str(i)+"========")
	Re_g = (G_g*D_h)/u_g
	Re_a = (G_a*D_h)/u_a
	jH_g,Cf_g = get_jH_and_Cf(Re_g)
	jH_a,Cf_a = get_jH_and_Cf(Re_a)

	h_g = jH_g*G_g*C_p_g*(Pr_g**(-2.0/3.0))
	h_a = jH_a*G_a*C_p_a*(Pr_a**(-2.0/3.0))

	m_g = np.sqrt( (h_g*wettedPerimeter)/(fin_k*A_k) )
	m_a = np.sqrt( (h_a*wettedPerimeter)/(fin_k*A_k) )
	n_f_g = np.tanh(m_g*finHeight)/(m_g*finHeight)
	n_f_a = np.tanh(m_a*finHeight)/(m_a*finHeight)
	n_o_g = 1.0-(1.0-n_f_g)*finAreaRatio
	n_o_a = 1.0-(1.0-n_f_a)*finAreaRatio

	R_f = finThickness/fin_k
	U = None
	if (i==0):
		U = 1.0/( 1.0/(n_o_g*h_g) + (R_f/n_o_g) + (alpha/alpha)*( R_f/n_o_a + 1.0/(n_o_a*h_a) ) )
	else:
		numPlates = (L-wallThickness)/(plateSpacing+wallThickness)
		A_w = L_g*L_a*numPlates
		U = 1.0/( 1.0/(n_o_g*h_g) + (R_f/n_o_g) + wallThickness/(wall_k*A_w) + (alpha/alpha)*( R_f/n_o_a + 1.0/(n_o_a*h_a) ) )

	print("U:"+str(U))
	A_g = (NTU*C_min)/U
	A_a = (alpha/alpha)*A_g #will be the same for this Hx
	A_o_g = m_dot_g/G_g #free flow area
	A_o_a = m_dot_a/G_a
	sigma = alpha*D_h/4.0
	A_fr_g = A_o_g/sigma
	A_fr_a = A_o_a/sigma

	L_g = L1 = (D_h*A_g)/(4.0*A_o_g)
	L_a = L2 = (D_h*A_a)/(4.0*A_o_a)
	L = L3 = A_fr_g/L2

	dP_g_calculated = ( (G_g**2)/(2.0*Rho_g_i) ) * ( 1.0-(sigma**2) + K_c + 2.0*(Rho_g_i/Rho_g_o-1.0) + Cf_g*L_g/r_h*Rho_g_i*(1.0/Rho_g_m) - (1.0-sigma**2-K_e)*(Rho_g_i/Rho_g_o) )
	dP_a_calculated = ( (G_a**2)/(2.0*Rho_a_i) ) * ( 1.0-(sigma**2) + K_c + 2.0*(Rho_a_i/Rho_a_o-1.0) + Cf_a*L_a/r_h*Rho_a_i*(1.0/Rho_a_m) - (1.0-sigma**2-K_e)*(Rho_a_i/Rho_a_o) )

	print("L_g:"+str(L_g))
	print("L_a:"+str(L_a))
	print("L:"+str(L))
	print("dP_g_calculated:"+str(dP_g_calculated))
	print("dP_a_calculated:"+str(dP_a_calculated))

	G_g = np.sqrt( (2.0*Rho_g_i*dP_g) / ( 1.0 - sigma**2 + K_c + 2.0*(Rho_g_i/Rho_g_o-1.0) + Cf_g*L_g/r_h*Rho_g_i*(1.0/Rho_g_m) - (1.0-sigma**2-K_e)*(Rho_g_i/Rho_g_o) ) )
	G_a = np.sqrt( (2.0*Rho_a_i*dP_a) / ( 1.0 - sigma**2 + K_c + 2.0*(Rho_a_i/Rho_a_o-1.0) + Cf_a*L_a/r_h*Rho_a_i*(1.0/Rho_a_m) - (1.0-sigma**2-K_e)*(Rho_a_i/Rho_a_o) ) )

	print(G_g)
	print(G_a)





