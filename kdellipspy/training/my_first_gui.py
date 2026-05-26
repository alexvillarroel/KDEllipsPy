import sys
import os
from PyQt6.QtWidgets import QApplication, QMainWindow
from PyQt6.QtGui import QIcon

# Rutas
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ICONO = os.path.join(BASE_DIR, "kdellipspy.ico")
SCRIPT = os.path.abspath(__file__)
DESKTOP_DIR = os.path.expanduser("~/.local/share/applications")
DESKTOP_FILE = os.path.join(DESKTOP_DIR, "kdellipspy.desktop")

def crear_desktop():
    os.makedirs(DESKTOP_DIR, exist_ok=True)
    contenido = f"""[Desktop Entry]
Name=KDEllipspy
Type=Application
Exec=python3 {SCRIPT}
Icon={ICONO}
"""
    with open(DESKTOP_FILE, "w") as f:
        f.write(contenido)
    os.system("update-desktop-database ~/.local/share/applications/")

# Solo lo crea si no existe o si el ícono cambió
if not os.path.exists(DESKTOP_FILE):
    crear_desktop()

app = QApplication(sys.argv)
app.setDesktopFileName("kdellipspy")

icono = QIcon(ICONO)
app.setWindowIcon(icono)

ventana = QMainWindow()
ventana.setWindowTitle("KDEllipspy")
ventana.setWindowIcon(icono)
ventana.resize(800, 600)
ventana.show()

sys.exit(app.exec())