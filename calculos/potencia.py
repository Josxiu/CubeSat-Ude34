# Generacion solar 1U, celdas en 4 caras laterales (+-X, +-Y); -Z camara (nadir), +Z antena
import numpy as np
S=1361; Acell=30.18e-4; ncell=2; eta=0.28   # 2 celdas TJ ~30 cm2 por cara (p.ej. Azur 3G30C), 28% conservador
k_loss=0.90*0.95*0.97  # MPPT/conv * temperatura * cableado/desajuste
Tmin=95.65; Tsun=61.0
Pface=S*Acell*ncell*eta*k_loss
print(f"Potencia por cara a incidencia normal: {Pface:.2f} W")
for beta in (15,22,30):
    b=np.radians(beta); u=np.linspace(0,2*np.pi,20000,endpoint=False)
    # marco LVLH: Z nadir, X velocidad, Y normal a la orbita. Sol en marco orbital:
    s_orb=np.array([np.cos(b),0,np.sin(b)])  # en marco (hacia nodo, ..., normal)
    # posicion angular u medida desde el punto subsolar proyectado
    r=np.vstack([np.cos(u),np.sin(u),0*u]).T     # radial
    v=np.vstack([-np.sin(u),np.cos(u),0*u]).T    # velocidad
    nrm=np.array([0,0,1.0])
    sx=v@s_orb; sy=np.full_like(u,s_orb@nrm); sz=-(r@s_orb)
    sunlit=np.degrees(np.abs(u))  # placeholder
    # eclipse: lado nocturno cuando r.s<0 y dentro del cilindro de sombra
    a=6928.137; Re=6378.137
    rs=r@s_orb; perp=np.sqrt(np.maximum(0,1-rs**2))*a
    ecl=(rs<0)&(perp<Re)
    P=Pface*(np.maximum(sx,0)+np.maximum(-sx,0)+np.maximum(sy,0)+np.maximum(-sy,0))
    P[ecl]=0
    Porb=P.mean(); E=Porb*Tmin/60
    print(f"beta {beta}: nadir -> P media orbital {Porb:.2f} W, E/orbita {E:.2f} Wh, P media en sol {P[~ecl].mean():.2f} W, eclipse {ecl.mean()*Tmin:.1f} min")
# modo tumbling/seguro: area proyectada media de 4 caras = 1.0 cara
Ptum=Pface*1.0*Tsun/Tmin
print(f"tumbling: P media orbital {Ptum:.2f} W, E/orbita {Ptum*Tmin/60:.2f} Wh")
