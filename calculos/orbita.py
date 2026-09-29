# Calculos de soporte Lab003 - orbita SSO 550 km, pasos sobre Medellin, eclipse y generacion
import numpy as np
mu=398600.4418; Re=6378.137; h=550.0; a=Re+h; J2=1.08263e-3
T=2*np.pi*np.sqrt(a**3/mu)
# inclinacion SSO
n=np.sqrt(mu/a**3); rate=2*np.pi/(365.2422*86400)
i=np.arccos(-rate/(1.5*n*J2*(Re/a)**2))
print(f"Periodo {T/60:.2f} min, orbitas/dia {86400/T:.2f}, inclinacion SSO {np.degrees(i):.2f} deg")
# geometria de visibilidad
rho=np.arcsin(Re/a)
for el in (5,10,15):
    e=np.radians(el); eta=np.arcsin(np.sin(rho)*np.cos(e)); lam=np.pi/2-e-eta
    d=Re*np.sin(lam)/np.sin(eta)
    print(f"elev min {el}: angulo central {np.degrees(lam):.1f} deg, rango max {d:.0f} km")
# eclipse vs beta (orbita circular, sombra cilindrica)
for beta in (0,15,22,30):
    b=np.radians(beta)
    val=np.sqrt(h*(h+2*Re))/(a*np.cos(b))
    fe=np.degrees(np.arccos(val))/180 if val<1 else 0
    print(f"beta {beta}: eclipse {fe*T/60:.1f} min, sol {(1-fe)*T/60:.1f} min")
# Simulacion de pasos (30 dias, paso 5 s), LTAN 10:30 (nodo descendente)
lat0,lon0=np.radians(6.25),np.radians(-75.57)
we=7.2921159e-5; dt=5.0; N=int(30*86400/dt); t=np.arange(N)*dt
# RAAN tal que nodo descendente cruza a 10:30 hora solar local en t=0 (aprox: sol en x inercial al inicio)
raan0=np.radians(-22.5+180)  # nodo ascendente a 22:30 LST
raan=raan0+rate*t
u=n*t
x=a*(np.cos(raan)*np.cos(u)-np.sin(raan)*np.sin(u)*np.cos(i))
y=a*(np.sin(raan)*np.cos(u)+np.cos(raan)*np.sin(u)*np.cos(i))
z=a*(np.sin(u)*np.sin(i))
th=we*t  # GMST ~ 0 al inicio (sol en +x, medianoche en Greenwich no importa para estadistica)
xe=np.cos(th)*x+np.sin(th)*y; ye=-np.sin(th)*x+np.cos(th)*y; ze=z
gs=Re*np.array([np.cos(lat0)*np.cos(lon0),np.cos(lat0)*np.sin(lon0),np.sin(lat0)])
r=np.vstack([xe,ye,ze]).T-gs
up=gs/np.linalg.norm(gs)
el=np.degrees(np.arcsin((r@up)/np.linalg.norm(r,axis=1)))
for emin in (10,15):
    vis=el>=emin
    edges=np.diff(vis.astype(int))
    starts=np.where(edges==1)[0]; ends=np.where(edges==-1)[0]
    if ends[0]<starts[0]: ends=ends[1:]
    m=min(len(starts),len(ends)); dur=(ends[:m]-starts[:m])*dt/60
    print(f"emin {emin}: pasos/dia {m/30:.2f}, duracion media {dur.mean():.1f} min, max {dur.max():.1f} min, contacto/dia {dur.sum()/30:.1f} min")
