import sys
from pathlib import Path
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QFileDialog, QTabWidget, 
                             QTextEdit, QTableWidget, QTableWidgetItem, QLabel, QHeaderView)
from PyQt6.QtCore import Qt

# Ensure the parent directory is in sys.path so we can import kdellipspy modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from kdellipspy.config_parser import ConfigParser

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("KDEllipsPy Inversion GUI")
        self.resize(900, 600)
        self.current_config_path = None
        self.cfg = None
        
        self._init_ui()
        
    def _init_ui(self):
        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)
        
        # Top control panel
        top_panel = QHBoxLayout()
        self.btn_load = QPushButton("Load input.ctl")
        self.btn_load.clicked.connect(self.load_input_file)
        self.lbl_status = QLabel("No file loaded.")
        top_panel.addWidget(self.btn_load)
        top_panel.addWidget(self.lbl_status)
        top_panel.addStretch()
        layout.addLayout(top_panel)
        
        # Tabs
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)
        
        # Tab 1: General Config
        self.tab_config = QWidget()
        self.tab_config_layout = QVBoxLayout(self.tab_config)
        self.txt_config_view = QTextEdit()
        self.txt_config_view.setReadOnly(True)
        self.tab_config_layout.addWidget(self.txt_config_view)
        self.tabs.addTab(self.tab_config, "Configuration")
        
        # Tab 2: Stations
        self.tab_stations = QWidget()
        self.tab_stations_layout = QVBoxLayout(self.tab_stations)
        self.table_stations = QTableWidget()
        self.table_stations.setColumnCount(7)
        self.table_stations.setHorizontalHeaderLabels(
            ["Name", "Lat", "Lon", "Elev", "Use N", "Use E", "Use Z"]
        )
        self.table_stations.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tab_stations_layout.addWidget(self.table_stations)
        self.tabs.addTab(self.tab_stations, "Stations")
        
        # Tab 3: Run / Console
        self.tab_run = QWidget()
        self.tab_run_layout = QVBoxLayout(self.tab_run)
        self.btn_run = QPushButton("Run Inversion")
        self.btn_run.setEnabled(False)
        self.txt_console = QTextEdit()
        self.txt_console.setReadOnly(True)
        self.tab_run_layout.addWidget(self.btn_run)
        self.tab_run_layout.addWidget(self.txt_console)
        self.tabs.addTab(self.tab_run, "Run & Log")

    def load_input_file(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self, "Select input.ctl", "", "Control Files (*.ctl *.*)"
        )
        if file_name:
            self.current_config_path = file_name
            self.lbl_status.setText(f"Loaded: {os.path.basename(file_name)}")
            try:
                self.cfg = ConfigParser(file_name)
                self.populate_ui_from_config()
                self.btn_run.setEnabled(True)
                self.txt_console.append(f"Successfully loaded config: {file_name}")
            except Exception as e:
                self.txt_console.append(f"Error loading config: {str(e)}")
                
    def populate_ui_from_config(self):
        if not self.cfg: return
        
        # Populate config text view (summary)
        config_summary = f"Event Name: {self.cfg.source_position.event_name}\n"
        config_summary += f"Fault Dimensions: {self.cfg.fault_plane.lx}x{self.cfg.fault_plane.ly} m\n"
        config_summary += f"Algorithm: {'MCMC' if self.cfg.inversion_process.algorithm_type == 1 else 'NA'}\n"
        
        tw = self.cfg.inversion_process.misfit_time_window
        tw_str = f"{tw} s" if tw > 0 else "Full Signal (0.0)"
        config_summary += f"Misfit Time Window: {tw_str}\n"
        
        self.txt_config_view.setText(config_summary)
        
        # Populate Stations table
        stations = self.cfg.stations.stations
        self.table_stations.setRowCount(len(stations))
        for i, st in enumerate(stations):
            self.table_stations.setItem(i, 0, QTableWidgetItem(st.name))
            self.table_stations.setItem(i, 1, QTableWidgetItem(str(st.latitude)))
            self.table_stations.setItem(i, 2, QTableWidgetItem(str(st.longitude)))
            self.table_stations.setItem(i, 3, QTableWidgetItem(str(st.height)))
            
            # Checkboxes for N, E, Z
            chk_n = QTableWidgetItem("N")
            chk_n.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            chk_n.setCheckState(Qt.CheckState.Checked if st.use_n else Qt.CheckState.Unchecked)
            self.table_stations.setItem(i, 4, chk_n)
            
            chk_e = QTableWidgetItem("E")
            chk_e.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            chk_e.setCheckState(Qt.CheckState.Checked if st.use_e else Qt.CheckState.Unchecked)
            self.table_stations.setItem(i, 5, chk_e)
            
            chk_z = QTableWidgetItem("Z")
            chk_z.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            chk_z.setCheckState(Qt.CheckState.Checked if st.use_z else Qt.CheckState.Unchecked)
            self.table_stations.setItem(i, 6, chk_z)

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
