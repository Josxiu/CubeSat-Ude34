# Oportunidades de imagen sobre la zona cafetera central (aprox. lat 1.5-7.0 N, lon -76.8 a -74.5)
import numpy as np
mu=398600.4418; Re=6378.137; a=Re+550; n=np.sqrt(mu/a**3); i=np.radians(97.59)
we=7.2921159e-5; rate=2*np.pi/(365.2422*86400)
dt=2.0; days=180; t=np.arange(0,days*86400,dt)
# Sol fijo en +X inercial (aprox. suficiente); nodo descendente a las 10:30 LST
raan=np.radians(180-22.5)+rate*t; u=n*t
x=np.cos(raan)*np.cos(u)-np.sin(raan)*np.sin(u)*np.cos(i)
y=np.sin(raan)*np.cos(u)+np.cos(raan)*np.sin(u)*np.cos(i)
z=np.sin(u)*np.sin(i)
sun_ang=rate*t  # el sol avanza en ascension recta
lst_sun=np.arctan2(y,x)-sun_ang  # angulo horario respecto al sol
day=np.cos(lst_sun)>0.3  # lado diurno (hora local ~7-17h)
th=we*t
lon=np.degrees(np.arctan2(-np.sin(th)*x+np.cos(th)*y, np.cos(th)*x+np.sin(th)*y))
lat=np.degrees(np.arcsin(z))
for halfswath_km in (30,52):
    m=halfswath_km/111.0
    inbox=(lat>1.5)&(lat<7.0)&(lon>-76.8-m)&(lon<-74.5+m)&day
    d=np.unique((t[inbox]//86400).astype(int))
    print(f"semi-franja {halfswath_km} km: dias con oportunidad {len(d)}/{days} -> 1 cada {days/len(d):.1f} dias")
