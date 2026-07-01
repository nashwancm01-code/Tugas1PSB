import streamlit as st
import math
import cmath
import matplotlib.pyplot as plt

# Konfigurasi halaman agar lebar (menyerupai layout Delphi)
st.set_page_config(page_title="Butterworth Filter Orde 2", layout="wide")

st.title("Desain Filter Digital Butterworth Orde 2")
st.write("Dihitung murni menggunakan *Bilinear Transformation* tanpa library SciPy/Numpy.")

# --- INPUT UTAMA ---
Ts = st.number_input("Tsampling (s):", value=0.01, step=0.001, format="%.4f")
fs = 1.0 / Ts
nyquist = fs / 2.0
st.info(f"Frekuensi Sampling (Fs) = {fs} Hz | Batas Nyquist = {nyquist} Hz")

st.markdown("---")

# --- FUNGSI MATEMATIKA MURNI (TANPA NUMPY) ---
def get_frequency_response(a0, a1, a2, b1, b2, fs_val):
    """Menghitung respon frekuensi H(f) murni menggunakan math dan cmath"""
    f_vals = []
    mag_vals = []
    nyq = fs_val / 2.0
    step = nyq / 100.0  # Resolusi 100 titik
    
    for i in range(101):
        f = i * step
        w = 2.0 * math.pi * f
        wTs = w / fs_val
        
        # Transformasi Z ke domain frekuensi: z = e^(j*w*Ts)
        # z^-1 = cos(-wTs) + j*sin(-wTs)
        z1 = complex(math.cos(-wTs), math.sin(-wTs))
        z2 = complex(math.cos(-2*wTs), math.sin(-2*wTs))
        
        # H(z) = (a0 + a1*z^-1 + a2*z^-2) / (1 + b1*z^-1 + b2*z^-2)
        num = a0 + a1*z1 + a2*z2
        den = 1.0 + b1*z1 + b2*z2
        
        H = 0 if abs(den) == 0 else num / den
        
        f_vals.append(f)
        mag_vals.append(abs(H))
        
    return f_vals, mag_vals

def plot_chart(f_vals, mag_vals, title):
    """Fungsi pembantu untuk memplot grafik menggunakan Matplotlib"""
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(f_vals, mag_vals, color='blue')
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Frekuensi [Hz]", fontsize=8)
    ax.set_ylabel("|H(f)|", fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.7)
    # Limit y-axis sedikit di atas 1 agar grafik rapi
    ax.set_ylim(0, 1.1)
    fig.tight_layout()
    return fig

# --- LAYOUT 2x2 SEPERTI DELPHI ---
col1, col2 = st.columns(2)
col3, col4 = st.columns(2)

# ==========================================
# 1. LOW PASS FILTER (LPF) - KIRI ATAS
# ==========================================
with col1:
    st.subheader("Low Pass Filter (LPF)")
    fc_lpf = st.slider("Fc LPF (Hz)", min_value=1.0, max_value=nyquist-1.0, value=40.0, step=1.0)
    
    # Perhitungan Bilinear untuk LPF
    omega_lpf = math.tan(math.pi * fc_lpf * Ts)
    omega2_lpf = omega_lpf ** 2
    C_lpf = 1.0 + math.sqrt(2)*omega_lpf + omega2_lpf
    
    a0_lpf = omega2_lpf / C_lpf
    a1_lpf = 2.0 * a0_lpf
    a2_lpf = a0_lpf
    b1_lpf = 2.0 * (omega2_lpf - 1.0) / C_lpf
    b2_lpf = (1.0 - math.sqrt(2)*omega_lpf + omega2_lpf) / C_lpf
    
    # Tampilan Koefisien
    c1, c2 = st.columns(2)
    c1.text_input("b1", value=f"{b1_lpf:.15f}", key="lpf_b1")
    c1.text_input("b2", value=f"{b2_lpf:.15f}", key="lpf_b2")
    c2.text_input("a0", value=f"{a0_lpf:.15f}", key="lpf_a0")
    c2.text_input("a1", value=f"{a1_lpf:.15f}", key="lpf_a1")
    c2.text_input("a2", value=f"{a2_lpf:.15f}", key="lpf_a2")
    
    st.latex(rf"Y(n) = {-b1_lpf:.4f}Y(n-1) - {b2_lpf:.4f}Y(n-2) + {a0_lpf:.4f}X(n) + {a1_lpf:.4f}X(n-1) + {a2_lpf:.4f}X(n-2)")
    
    f_lpf, mag_lpf = get_frequency_response(a0_lpf, a1_lpf, a2_lpf, b1_lpf, b2_lpf, fs)
    st.pyplot(plot_chart(f_lpf, mag_lpf, "Respon Filter LPF"))

# ==========================================
# 2. HIGH PASS FILTER (HPF) - KANAN ATAS
# ==========================================
with col2:
    st.subheader("High Pass Filter (HPF)")
    fc_hpf = st.slider("Fc HPF (Hz)", min_value=1.0, max_value=nyquist-1.0, value=4.0, step=1.0)
    
    # Perhitungan Bilinear untuk HPF
    omega_hpf = math.tan(math.pi * fc_hpf * Ts)
    omega2_hpf = omega_hpf ** 2
    C_hpf = 1.0 + math.sqrt(2)*omega_hpf + omega2_hpf
    
    a0_hpf = 1.0 / C_hpf
    a1_hpf = -2.0 / C_hpf
    a2_hpf = 1.0 / C_hpf
    b1_hpf = 2.0 * (omega2_hpf - 1.0) / C_hpf
    b2_hpf = (1.0 - math.sqrt(2)*omega_hpf + omega2_hpf) / C_hpf
    
    # Tampilan Koefisien
    c1, c2 = st.columns(2)
    c1.text_input("b1", value=f"{b1_hpf:.15f}", key="hpf_b1")
    c1.text_input("b2", value=f"{b2_hpf:.15f}", key="hpf_b2")
    c2.text_input("a0", value=f"{a0_hpf:.15f}", key="hpf_a0")
    c2.text_input("a1", value=f"{a1_hpf:.15f}", key="hpf_a1")
    c2.text_input("a2", value=f"{a2_hpf:.15f}", key="hpf_a2")
    
    st.latex(rf"Y(n) = {-b1_hpf:.4f}Y(n-1) - {b2_hpf:.4f}Y(n-2) + {a0_hpf:.4f}X(n) + {a1_hpf:.4f}X(n-1) + {a2_hpf:.4f}X(n-2)")

    f_hpf, mag_hpf = get_frequency_response(a0_hpf, a1_hpf, a2_hpf, b1_hpf, b2_hpf, fs)
    st.pyplot(plot_chart(f_hpf, mag_hpf, "Respon Filter HPF"))

st.markdown("---")

# ==========================================
# 3. BAND PASS FILTER (BPF) - KIRI BAWAH
# ==========================================
with col3:
    st.subheader("Band Pass Filter (BPF)")
    
    # Slider dibuat terpisah
    fc_low_bpf = st.slider("Fc Low BPF (Hz)", min_value=1.0, max_value=nyquist-2.0, value=10.0, step=1.0)
    fc_high_bpf = st.slider("Fc High BPF (Hz)", min_value=2.0, max_value=nyquist-1.0, value=25.0, step=1.0)
    
    # Validasi logika frekuensi
    if fc_low_bpf >= fc_high_bpf:
        st.warning("⚠️ Nilai **Fc High** harus lebih besar dari **Fc Low**!")
    else:
        # Perhitungan Bilinear untuk BPF
        O_L_bpf = math.tan(math.pi * fc_low_bpf * Ts)
        O_H_bpf = math.tan(math.pi * fc_high_bpf * Ts)
        BW_bpf = O_H_bpf - O_L_bpf
        w02_bpf = O_L_bpf * O_H_bpf
        C_bpf = 1.0 + BW_bpf + w02_bpf
        
        a0_bpf = BW_bpf / C_bpf
        a1_bpf = 0.0
        a2_bpf = -BW_bpf / C_bpf
        b1_bpf = 2.0 * (w02_bpf - 1.0) / C_bpf
        b2_bpf = (1.0 - BW_bpf + w02_bpf) / C_bpf
        
        # Tampilan Koefisien
        c1, c2 = st.columns(2)
        c1.text_input("b1", value=f"{b1_bpf:.15f}", key="bpf_b1")
        c1.text_input("b2", value=f"{b2_bpf:.15f}", key="bpf_b2")
        c2.text_input("a0", value=f"{a0_bpf:.15f}", key="bpf_a0")
        c2.text_input("a1", value=f"{a1_bpf:.15f}", key="bpf_a1")
        c2.text_input("a2", value=f"{a2_bpf:.15f}", key="bpf_a2")
        
        f_bpf, mag_bpf = get_frequency_response(a0_bpf, a1_bpf, a2_bpf, b1_bpf, b2_bpf, fs)
        st.pyplot(plot_chart(f_bpf, mag_bpf, "Respon Filter BPF"))

# ==========================================
# 4. BAND STOP FILTER (BSF) - KANAN BAWAH
# ==========================================
with col4:
    st.subheader("Band Stop Filter (BSF)")
    
    # Slider dibuat terpisah
    fc_low_bsf = st.slider("Fc Low BSF (Hz)", min_value=1.0, max_value=nyquist-2.0, value=10.0, step=1.0)
    fc_high_bsf = st.slider("Fc High BSF (Hz)", min_value=2.0, max_value=nyquist-1.0, value=35.0, step=1.0)
    
    # Validasi logika frekuensi
    if fc_low_bsf >= fc_high_bsf:
        st.warning("⚠️ Nilai **Fc High** harus lebih besar dari **Fc Low**!")
    else:
        # Perhitungan Bilinear untuk BSF
        O_L_bsf = math.tan(math.pi * fc_low_bsf * Ts)
        O_H_bsf = math.tan(math.pi * fc_high_bsf * Ts)
        BW_bsf = O_H_bsf - O_L_bsf
        w02_bsf = O_L_bsf * O_H_bsf
        C_bsf = 1.0 + BW_bsf + w02_bsf
        
        a0_bsf = (1.0 + w02_bsf) / C_bsf
        a1_bsf = 2.0 * (w02_bsf - 1.0) / C_bsf
        a2_bsf = (1.0 + w02_bsf) / C_bsf
        b1_bsf = 2.0 * (w02_bsf - 1.0) / C_bsf
        b2_bsf = (1.0 - BW_bsf + w02_bsf) / C_bsf
        
        # Tampilan Koefisien
        c1, c2 = st.columns(2)
        c1.text_input("b1", value=f"{b1_bsf:.15f}", key="bsf_b1")
        c1.text_input("b2", value=f"{b2_bsf:.15f}", key="bsf_b2")
        c2.text_input("a0", value=f"{a0_bsf:.15f}", key="bsf_a0")
        c2.text_input("a1", value=f"{a1_bsf:.15f}", key="bsf_a1")
        c2.text_input("a2", value=f"{a2_bsf:.15f}", key="bsf_a2")
        
        f_bsf, mag_bsf = get_frequency_response(a0_bsf, a1_bsf, a2_bsf, b1_bsf, b2_bsf, fs)
        st.pyplot(plot_chart(f_bsf, mag_bsf, "Respon Filter BSF"))
