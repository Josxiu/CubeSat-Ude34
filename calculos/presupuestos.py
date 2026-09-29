# Presupuestos preliminares Lab003 - CubeSat UD34
# Enlace UHF, carga util (GSD/datos), almacenamiento, energia, masa y volumen.
# Todas las cifras del informe v2 salen de este script (python presupuestos.py).
import numpy as np

# ---------------- Orbita y contactos (ver orbita.py / cobertura.py) ----------------
T_orb_min = 95.65          # periodo orbital [min]
orb_dia = 1440 / T_orb_min
t_ecl_min = 35.2           # eclipse peor caso (beta 15 deg)
contacto_min_dia = 15.6    # contacto total/dia sobre Medellin, elev >= 10 deg (simulado)
pasos_dia = 2.37
paso_max_min = 7.9
dias_entre_oport = 5.5     # oportunidades de imagen zona cafetera central (semi-franja 52 km)

# ---------------- Carga util ----------------
h = 550e3
pix = 2.2e-6               # pitch del sensor (clase 5 MP, 2592 x 1944)
nx, ny = 2592, 1944
GSD_obj = 40.0
f = h * pix / GSD_obj
swath = nx * GSD_obj / 1e3
v_suelo = 7.586 * 6378.137 / 6928.137   # km/s velocidad de la traza
t_exp_max = 0.5 * GSD_obj / (v_suelo * 1e3)   # arrastre <= 0.5 px
bits = 10
bandas = 4
fila_banda = ny // bandas
L_franja = fila_banda * GSD_obj / 1e3          # km por franja de filtro
dt_frame = L_franja / v_suelo                  # s entre cuadros (pushframe)
frame_MB = nx * ny * bits / 8 / 1e6
L_tira_km = 470                                 # recorrido N-S de la zona cafetera central (~1.5-5.8 N)
n_frames = int(np.ceil((L_tira_km + ny * GSD_obj / 1e3) / L_franja))
MB_oport = n_frames * frame_MB
ang_apunt = np.degrees(np.arctan((swath / 2 - 10) / 550))
print("=== CARGA UTIL ===")
print(f"f = {f*1e3:.1f} mm, franja = {swath:.1f} km, v_suelo = {v_suelo:.2f} km/s")
print(f"t_exp max (0.5 px) = {t_exp_max*1e3:.2f} ms; franja por banda = {L_franja:.1f} km; cuadro cada {dt_frame:.2f} s")
print(f"cuadro crudo = {frame_MB:.2f} MB; cuadros por oportunidad = {n_frames}; datos crudos/oportunidad = {MB_oport:.0f} MB")
print(f"apuntamiento requerido para dejar el objetivo a >=10 km del borde: {ang_apunt:.1f} deg")
# difraccion
D = f / 2.8
airy = 1.22 * 0.84e-6 / D * f
print(f"f/2.8 -> D = {D*1e3:.1f} mm; radio de Airy a 840 nm = {airy*1e6:.2f} um ({airy/pix:.2f} px)")

# ---------------- Enlace UHF (435 MHz) ----------------
k = -228.6
fc = 436e6
lam = 3e8 / fc
def fspl(d_km):
    return 20 * np.log10(4 * np.pi * d_km * 1e3 / lam)
print("\n=== ENLACE UHF DESCENDENTE ===")
Ptx_dBW = 0.0          # 30 dBm (1 W)
G_sat = 0.0            # dipolo/turnstile, peor direccion util
L_linea = 1.0
L_pol = 3.0            # lineal (sat) -> circular (tierra)
L_atm_iono = 2.0       # atmosfera + centelleo ionosferico (latitud ecuatorial)
L_apunt = 0.5
G_gs = 14.0            # Yagi cruzada ~14 dBic
Tsys = 1000.0          # K, ruido urbano UHF conservador
GT = G_gs - 10 * np.log10(Tsys)
EbN0_req = 10.0        # GFSK/GMSK + FEC (RS+conv), BER 1e-5, con perdida de implementacion
for el, d in ((10, 1816), (30, 1002), (90, 550)):
    CN0 = Ptx_dBW + G_sat - L_linea - fspl(d) - L_pol - L_atm_iono - L_apunt + GT - k
    s = f"elev {el:2d} deg (d={d} km): FSPL={fspl(d):.1f} dB, C/N0={CN0:.1f} dBHz"
    for R in (9600, 19200, 38400):
        s += f" | {R/1000:.1f}k margen {CN0 - 10*np.log10(R) - EbN0_req:+.1f} dB"
    print(s)
print("\n=== ENLACE UHF ASCENDENTE ===")
P_gs = 10 * np.log10(50)   # 50 W
Prx = P_gs + 30 + G_gs - 2 - fspl(1816) - L_pol - L_atm_iono - L_linea + G_sat
print(f"Prx en el satelite (10 deg) = {Prx:.1f} dBm vs sensibilidad AX100 -122 dBm @ 5 kbaud -> margen {Prx+122:.1f} dB")

# ---------------- Datos ----------------
print("\n=== DATOS ===")
eff_proto = 0.75       # overhead de trama/FEC/reintentos
disp = 0.8             # fraccion de pasos realmente usados
for R in (9600, 19200):
    MBd = R * contacto_min_dia * 60 * eff_proto * disp / 8 / 1e6
    print(f"{R/1000:.1f} kbps: capacidad util {MBd:.2f} MB/dia, {MBd*dias_entre_oport:.1f} MB por ciclo de {dias_entre_oport} dias")
comp = 2.0             # compresion sin perdidas (CCSDS 123 / JPEG-LS), tipica 2:1
roi_km = 50
roi_MB = (roi_km * 1e3 / GSD_obj) ** 2 * bandas * bits / 8 / 1e6 / comp
quick_MB = MB_oport / 256 / comp  # quicklook decimado 16x16
tm_MB_dia = 64 * 1440 / 1e6       # 64 B/min de telemetria de mision
print(f"ROI {roi_km}x{roi_km} km, 4 bandas, comprimida: {roi_MB:.2f} MB; quicklook por oportunidad {quick_MB:.2f} MB; TM {tm_MB_dia:.3f} MB/dia")
# almacenamiento
n_guardadas = 10
alm = n_guardadas * MB_oport / 1e3
print(f"Almacenamiento crudo de {n_guardadas} oportunidades = {alm:.2f} GB -> con margen 50% = {alm/0.5:.2f} GB")
buf = 2 * nx * ny * 2 / 1e6     # doble buffer, 16 bit por pixel
print(f"Buffer RAM doble cuadro (16 bit) = {buf:.1f} MB -> SDRAM >= {buf/0.5:.0f} MB (margen 50%)")
# bus interno
thr_req = 19.2e3
can_eff = 0.45e6
print(f"CAN 1 Mbps efectivo ~{can_eff/1e3:.0f} kbps vs {thr_req/1e3:.1f} kbps requeridos -> uso {thr_req/can_eff*100:.1f}%")

# ---------------- Energia ----------------
print("\n=== ENERGIA (por orbita, Wh) ===")
To = T_orb_min / 60
eta_conv = 0.88
tx_min_orb = contacto_min_dia * 0.6 / orb_dia      # Tx ~60% del tiempo de paso
cargas = [
    # nombre, P_nom [W], P_max [W], ciclo nominal (0-1), fuente
    ("OBC (MCU Cortex-M4/M33 clase baja potencia + FRAM + NOR)", 0.10, 0.20, 1.0),
    ("UHF Rx (clase AX100: 55 mA @3.3 V)", 0.18, 0.40, 1.0),
    ("UHF Tx 1 W RF (800 mA @3.3 V)", 2.64, 3.30, tx_min_orb / T_orb_min),
    ("ADCS sensores (mag + giro + 6 sol grueso)", 0.10, 0.15, 1.0),
    ("ADCS magnetorquers (3 ejes)", 0.20, 0.60, 0.5),
    ("ADCS rueda de momento (eje pitch)", 0.30, 1.00, 1.0),
    ("Carga util: sensor + procesador + SDRAM + eMMC", 1.70, 2.50, (10 / (dias_entre_oport * orb_dia)) / T_orb_min),
    ("Calefactor de bateria (nominal 0; peor caso ver abajo)", 1.00, 1.00, 0.0),
]
Eload = 0.0
for n, P, Pm, dc in cargas:
    E = P * dc * To
    Eload += E
    print(f"{n:58s} Pnom {P:4.2f} Pmax {Pm:4.2f} ciclo {dc*100:6.2f}% E {E:.3f}")
E_eps = 0.10 * To      # autoconsumo EPS (MPPT, monitores)
E_tot = Eload / eta_conv + E_eps
print(f"Cargas {Eload:.2f} Wh -> con conversion ({eta_conv:.0%}) + EPS {E_eps:.2f} = {E_tot:.2f} Wh/orbita")
for beta, Eg in ((15, 1.81), (22, 2.00), (30, 2.20)):
    print(f"  beta {beta}: generacion {Eg:.2f} Wh -> margen {(Eg-E_tot)/Eg*100:+.0f}%")
# peor orbita: paso completo + imagen + calefactor 20% del eclipse
E_peor = E_tot + (2.64 * paso_max_min / 60 + 1.70 * 10 / 60 + 1.0 * 0.2 * t_ecl_min / 60) / eta_conv
print(f"Peor orbita (paso 7.9 min Tx + 10 min carga util + calefactor): {E_peor:.2f} Wh -> deficit {E_peor-1.81:.2f} Wh cubierto por bateria")
E_ecl = (E_tot / To) * t_ecl_min / 60
Ebat = 2 * 3.6 * 2.6 * 1.0   # 2S1P 18650 2600 mAh
print(f"Energia en eclipse nominal {E_ecl:.2f} Wh -> DoD {E_ecl/Ebat*100:.1f}% de {Ebat:.1f} Wh; peor orbita DoD {(E_ecl + (E_peor-E_tot))/Ebat*100:.1f}%")
# modo seguro
P_seg = (0.10 + 0.18 + 0.10 + 0.20 * 0.5) / eta_conv + 0.10
print(f"Modo seguro (OBC+Rx+sensores+B-dot, sin rueda ni carga util): {P_seg:.2f} W medio -> {P_seg*To:.2f} Wh/orbita vs 1.94 Wh tumbling")
print(f"Autonomia solo bateria en modo seguro (80% DoD): {0.8*Ebat/P_seg:.0f} h")
# comparacion con EnduroSat UHF II
E_endu = 1.2 * To
print(f"Comparacion: EnduroSat UHF II en Rx (1.2 W) = {E_endu:.2f} Wh/orbita, solo escuchando")

# ---------------- Masa ----------------
print("\n=== MASA (g) ===")
masa = [
    ("Estructura 1U", 110), ("Tornilleria, soportes y arnes", 70),
    ("4 paneles solares laterales (2 celdas TJ c/u)", 160),
    ("PCB-1 EPS", 60), ("Bateria 2x18650 + soporte", 115),
    ("PCB-2 OBC/ADCS/COM", 60), ("Modulo transceptor UHF", 25),
    ("Antena UHF desplegable (+Z)", 85), ("Rueda de momento", 60),
    ("3 magnetorquers", 23), ("Sensores de sol gruesos", 10),
    ("PCB-3 carga util (sensor + procesador)", 60), ("Optica f~30 mm + filtro de 4 bandas + montura", 120),
]
m = sum(x[1] for x in masa)
for n, g in masa:
    print(f"{n:48s} {g:5d}")
print(f"Total {m} g; con 20% de margen {m*1.2:.0f} g ({m*1.2/2000*100:.0f}% de 2.00 kg)")

# ---------------- Volumen (alturas del stack a lo largo de Z) ----------------
print("\n=== VOLUMEN (mm en Z) ===")
stack = [("Antena UHF desplegable (+Z), altura a confirmar con datasheet", 12), ("PCB-1 EPS + baterias 18650 horizontales", 25),
         ("PCB-2 OBC/ADCS/COM (transceptor montado como hija)", 15),
         ("PCB-3 carga util (sensor en cara -Z, procesador/SDRAM/eMMC atras)", 12),
         ("Optica f~30 mm (barril sobre el sensor)", 30)]
hz = sum(x[1] for x in stack)
for n, z in stack:
    print(f"{n:52s} {z:3d}")
print(f"Total {hz} mm de ~98 mm internos utiles (estimado) -> holgura {98-hz} mm")
