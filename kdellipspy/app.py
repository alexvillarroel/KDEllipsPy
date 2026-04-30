import sys
from PySide6.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget, 
                             QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, 
                             QFormLayout, QTextEdit, QFileDialog)
from PySide6.QtCore import Qt, QThread, Signal

# Importamos tus clases existentes
from config_parser import ConfigParser
from inversion_na import NAInversionModel, NAConfig

class InversionWorker(QThread):
    """Hilo para ejecutar la inversión sin bloquear la GUI."""
    finished = Signal(object)
    log_signal = Signal(str)

    def __init__(self, model, config):
        super().__init__()
        self.model = model
        self.config = config

    def run(self):
        # Aquí llamarías a run_na_search de tu clase NAInversionModel
        result = self.model.run_na_search(self.config)
        self.finished.emit(result)

class EventTab(QWidget):
    """Pestaña para los parámetros de la sección 1 y 2 del input.ctl."""
    def __init__(self, cfg):
        super().__init__()
        layout = QFormLayout()
        self.event_name = QLineEdit(cfg.source_position.event_name)
        self.lat = QLineEdit(str(cfg.source_position.latitude))
        self.lon = QLineEdit(str(cfg.source_position.longitude))
        
        layout.addRow("Nombre del Evento:", self.event_name)
        layout.addRow("Latitud:", self.lat)
        layout.addRow("Longitud:", self.lon)
        self.setLayout(layout)

class KDEllipsPyGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("KDEllipsPy - Kinematic Inversion")
        self.resize(1000, 700)
        self.cfg = None

        # Widget central y pestañas
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        # Consola de salida
        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setPlaceholderText("Logs de la inversión...")

        self.init_menu()

    def init_menu(self):
        menu = self.menuBar().addMenu("Archivo")
        load_action = menu.addAction("Cargar input.ctl")
        load_action.triggered.connect(self.load_project)

    def load_project(self):
        path, _ = QFileDialog.getOpenFileName(self, "Abrir Control File", "", "CTL Files (*.ctl)")
        if path:
            self.cfg = ConfigParser(path)
            self.setup_tabs()

    def setup_tabs(self):
        self.tabs.clear()
        self.tabs.addTab(EventTab(self.cfg), "🌍 Evento")
        # Aquí añadirías pestañas para StationsTab, VelocityTab, etc.
        
        # Pestaña de ejecución
        run_widget = QWidget()
        run_layout = QVBoxLayout()
        run_btn = QPushButton("🚀 Iniciar Inversión NA")
        run_btn.clicked.connect(self.start_inversion)
        run_layout.addWidget(run_btn)
        run_layout.addWidget(self.console)
        run_widget.setLayout(run_layout)
        self.tabs.addTab(run_widget, "⚙️ Ejecución")

    def start_inversion(self):
        if not self.cfg: return
        self.console.append("Iniciando búsqueda NA...")
        
        # Configuración de inversión desde tu backend
        inversion_model = NAInversionModel(self.cfg.filepath) 
        na_cfg = NAConfig(n_iterations=self.cfg.inversion_process.num_iterations)
        
        self.worker = InversionWorker(inversion_model, na_cfg)
        self.worker.finished.connect(self.on_finished)
        self.worker.start()

    def on_finished(self, result):
        self.console.append(f"Inversión terminada. Mejor misfit: {result.best_model.misfit:.4e}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = KDEllipsPyGUI()
    window.show()
    sys.exit(app.exec())