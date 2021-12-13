# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pyqtgraph as pg
import pandas as pd
import pprint as pp

import BFO.BFOator.Bfield_001_type1 as type1
import BFO.BFOator.Bfield_001_type2 as type2
import BFO.BFOator.Bfield_111_type1 as t1_111
import BFO.BFOator.Bfield_111_type1_P_ip as t1_111_ip

from qtpy import QtWidgets
from qtpy import uic
from scipy.optimize import curve_fit
from scipy.signal import argrelmin

def project_B(Bx, By, Bz, theta, phi):
    return np.sin(theta)*np.cos(phi)*Bx + np.sin(theta)*np.sin(phi)*By + np.cos(theta)*Bz 

def sinus(x, a, b, c):
    """
    Sine function for fitting.
    a : amplitude
    b: period
    c: phase shift (in rad)
    """
    return a*np.sin(2*np.pi*x/b+c)

                  
class BFOatorMainWindow(QtWidgets.QMainWindow):
    """ The main window for the BFO-ator.
    """
    
    def __init__(self):
        # Get the path to the *.ui file
        this_dir = os.path.dirname(__file__)
        ui_file = os.path.join(this_dir, 'ui_gui.ui')

        # Load it
        super(BFOatorMainWindow, self).__init__()
        uic.loadUi(ui_file, self)
        self.show()


class BFOator():

    def __init__(self):
        super().__init__()
        self.activate_gui()

        
    def show(self):
        self._mw.show()

        
    def activate_gui(self):
        """
        Sets up and connects everything.
        """
        
        self._mw = BFOatorMainWindow()

        self.set_up_tab_001t1()
        self.set_up_tab_001t2()
        self.set_up_tab_111()
        self.set_up_tab_111_ip()
        self.cycloid_type = ""


    def set_up_tab_001t1(self):
        """
        Activation of the 001, type 1 tab.
        """
        
        self._mw.k_001t1_comboBox.addItems(["k1 [1-10]", "k2 [-101]", "k3 [0-11]"])

        # connect signals
        self._mw.compute_001t1_pushButton.clicked.connect(self.set_up_and_calc_001t1)
        self._mw.help_001t1_pushButton.clicked.connect(self.display_help_001t1)
        self._mw.save_001t1_pushButton.clicked.connect(self.save_routine)
        self._mw.Bx_001t1_checkBox.stateChanged.connect(self.update_plot)
        self._mw.By_001t1_checkBox.stateChanged.connect(self.update_plot)
        self._mw.Bz_001t1_checkBox.stateChanged.connect(self.update_plot)
        self._mw.BNV_001t1_checkBox.stateChanged.connect(self.update_plot)

        # disable save when no data computed yet
        self._mw.save_001t1_pushButton.setEnabled(False)

        # set up plot
        self.xmin = self._mw.xmin_001t1_doubleSpinBox.value()
        self.xmax = self._mw.xmax_001t1_doubleSpinBox.value()
        self.nb_pts = self._mw.nb_pts_001t1_spinBox.value()

        self.r_array = np.linspace(self.xmin, self.xmax, self.nb_pts)*1e-9
        self.data = pd.DataFrame({"Bx": np.zeros(self.nb_pts),
                                  "By": np.zeros(self.nb_pts),
                                  "Bz": np.zeros(self.nb_pts),
                                  "BNV": np.zeros(self.nb_pts)},
                                 index=self.r_array)
        self.Bx_plot_001t1 = pg.PlotDataItem(self.r_array, self.data["Bx"].values,
                                             pen=pg.mkPen(color=(239, 83, 80), width=2))
        self.By_plot_001t1 = pg.PlotDataItem(self.r_array, self.data["By"].values,
                                             pen=pg.mkPen(color=(70, 25, 136), width=2))
        self.Bz_plot_001t1 = pg.PlotDataItem(self.r_array, self.data["Bz"].values,
                                             pen=pg.mkPen(color=(22, 173, 170), width=2))
        self.BNV_plot_001t1 = pg.PlotDataItem(self.r_array, self.data["BNV"].values,
                                              pen=pg.mkPen(color=(255, 165, 0), width=2))

        self._mw.type1_001_graphicsView.addItem(self.Bx_plot_001t1)
        self._mw.type1_001_graphicsView.addItem(self.By_plot_001t1)
        self._mw.type1_001_graphicsView.addItem(self.Bz_plot_001t1)
        self._mw.type1_001_graphicsView.addItem(self.BNV_plot_001t1)
        self._mw.type1_001_graphicsView.setLabel("bottom", "r", units="m")
        self._mw.type1_001_graphicsView.setLabel("left", "B", units="T")
        self._mw.type1_001_graphicsView.setBackground("w")
        legend = pg.LegendItem(colCount=4)
        legend.addItem(self.Bx_plot_001t1, name="  Bx")
        legend.addItem(self.By_plot_001t1, name="  By")
        legend.addItem(self.Bz_plot_001t1, name="  Bz")
        legend.addItem(self.BNV_plot_001t1, name="  BNV")
        legend.setParentItem(self._mw.type1_001_graphicsView.getPlotItem())
        legend.setPos(100, -10)
        

    def set_up_tab_001t2(self):
        """
        Activation of the 001, type 2 tab.
        """
        
        self._mw.k_001t2_comboBox.addItems(["k'1 [-211]", "k'2 [1-21]", "k'3 [11-2]"])

        # connect signals
        self._mw.compute_001t2_pushButton.clicked.connect(self.set_up_and_calc_001t2)
        self._mw.help_001t2_pushButton.clicked.connect(self.display_help_001t2)
        self._mw.save_001t2_pushButton.clicked.connect(self.save_routine)
        self._mw.k_001t2_comboBox.currentIndexChanged.connect(self.enable_disable_alpha)
        self._mw.Bx_001t2_checkBox.stateChanged.connect(self.update_plot)
        self._mw.By_001t2_checkBox.stateChanged.connect(self.update_plot)
        self._mw.Bz_001t2_checkBox.stateChanged.connect(self.update_plot)
        self._mw.BNV_001t2_checkBox.stateChanged.connect(self.update_plot)

        # disable save when no data computed yet
        self._mw.save_001t2_pushButton.setEnabled(False)
        self._mw.alpha_doubleSpinBox.setEnabled(False)

        # set up plot
        self.xmin = self._mw.xmin_001t2_doubleSpinBox.value()
        self.xmax = self._mw.xmax_001t2_doubleSpinBox.value()
        self.nb_pts = self._mw.nb_pts_001t2_spinBox.value()

        self.r_array = np.linspace(self.xmin, self.xmax, self.nb_pts)*1e-9
        self.data = pd.DataFrame({"Bx": np.zeros(self.nb_pts),
                                  "By": np.zeros(self.nb_pts),
                                  "Bz": np.zeros(self.nb_pts),
                                  "BNV": np.zeros(self.nb_pts)},
                                 index=self.r_array)
        self.Bx_plot_001t2 = pg.PlotDataItem(self.r_array, self.data["Bx"].values,
                                             pen=pg.mkPen(color=(239, 83, 80), width=2))
        self.By_plot_001t2 = pg.PlotDataItem(self.r_array, self.data["By"].values,
                                             pen=pg.mkPen(color=(70, 25, 136), width=2))
        self.Bz_plot_001t2 = pg.PlotDataItem(self.r_array, self.data["Bz"].values,
                                             pen=pg.mkPen(color=(22, 173, 170), width=2))
        self.BNV_plot_001t2 = pg.PlotDataItem(self.r_array, self.data["BNV"].values,
                                              pen=pg.mkPen(color=(255, 165, 0), width=2))
                                             

        self._mw.type2_001_graphicsView.addItem(self.Bx_plot_001t2)
        self._mw.type2_001_graphicsView.addItem(self.By_plot_001t2)
        self._mw.type2_001_graphicsView.addItem(self.Bz_plot_001t2)
        self._mw.type2_001_graphicsView.addItem(self.BNV_plot_001t2)
        self._mw.type2_001_graphicsView.setLabel("bottom", "r", units="m")
        self._mw.type2_001_graphicsView.setLabel("left", "B", units="T")
        self._mw.type2_001_graphicsView.setBackground("w")
        legend = pg.LegendItem(colCount=4)
        legend.addItem(self.Bx_plot_001t2, name="  Bx")
        legend.addItem(self.By_plot_001t2, name="  By")
        legend.addItem(self.Bz_plot_001t2, name="  Bz")
        legend.addItem(self.BNV_plot_001t2, name="  BNV")
        legend.setParentItem(self._mw.type2_001_graphicsView.getPlotItem())
        legend.setPos(100, -10)
        return


    def enable_disable_alpha(self):
        if self._mw.k_001t2_comboBox.currentText()=="k'3 [11-2]":
            self._mw.alpha_doubleSpinBox.setEnabled(True)
        else:
            self._mw.alpha_doubleSpinBox.setEnabled(False)


    def set_up_tab_111(self):
        """
        Activation of the 111, type 1 tab.
        """

        # connect signals
        self._mw.compute_111_pushButton.clicked.connect(self.set_up_and_calc_111)
        self._mw.help_111_pushButton.clicked.connect(self.display_help_111)
        self._mw.save_111_pushButton.clicked.connect(self.save_routine)
        self._mw.Bx_111_checkBox.stateChanged.connect(self.update_plot)
        self._mw.By_111_checkBox.stateChanged.connect(self.update_plot)
        self._mw.Bz_111_checkBox.stateChanged.connect(self.update_plot)
        self._mw.BNV_111_checkBox.stateChanged.connect(self.update_plot)

        # disable save when no data computed yet
        self._mw.save_111_pushButton.setEnabled(False)

        # set up plot
        self.xmin = self._mw.xmin_111_doubleSpinBox.value()
        self.xmax = self._mw.xmax_111_doubleSpinBox.value()
        self.nb_pts = self._mw.nb_pts_111_spinBox.value()

        self.r_array = np.linspace(self.xmin, self.xmax, self.nb_pts)*1e-9
        self.data = pd.DataFrame({"Bx": np.zeros(self.nb_pts),
                                  "By": np.zeros(self.nb_pts),
                                  "Bz": np.zeros(self.nb_pts),
                                  "BNV": np.zeros(self.nb_pts)},
                                 index=self.r_array)
        self.Bx_plot_111 = pg.PlotDataItem(self.r_array, self.data["Bx"].values,
                                           pen=pg.mkPen(color=(239, 83, 80), width=2))
        self.By_plot_111 = pg.PlotDataItem(self.r_array, self.data["By"].values,
                                           pen=pg.mkPen(color=(70, 25, 136), width=2))         
        self.Bz_plot_111 = pg.PlotDataItem(self.r_array, self.data["Bz"].values,
                                           pen=pg.mkPen(color=(22, 173, 170), width=2))
        self.BNV_plot_111 = pg.PlotDataItem(self.r_array, self.data["BNV"].values,
                                            pen=pg.mkPen(color=(255, 165, 0), width=2))
                                             

        self._mw.type1_111_graphicsView.addItem(self.Bx_plot_111)
        self._mw.type1_111_graphicsView.addItem(self.By_plot_111)
        self._mw.type1_111_graphicsView.addItem(self.Bz_plot_111)
        self._mw.type1_111_graphicsView.addItem(self.BNV_plot_111)
        self._mw.type1_111_graphicsView.setLabel("bottom", "r", units="m")
        self._mw.type1_111_graphicsView.setLabel("left", "B", units="T")
        self._mw.type1_111_graphicsView.setBackground("w")

        legend = pg.LegendItem(colCount=4)
        legend.addItem(self.Bx_plot_111, name="  Bx")
        legend.addItem(self.By_plot_111, name="  By")
        legend.addItem(self.Bz_plot_111, name="  Bz")
        legend.addItem(self.BNV_plot_111, name="  BNV")
        legend.setParentItem(self._mw.type1_111_graphicsView.getPlotItem())
        legend.setPos(100, -10)
        
        return


    def set_up_tab_111_ip(self):
        """
        Activation of the 111, type 1 P in plane tab.
        """
        self._mw.k_111_ip_comboBox.addItems(["k1", "k2", "k3"])

        # connect signals
        self._mw.compute_111_ip_pushButton.clicked.connect(self.set_up_and_calc_111_ip)
        self._mw.help_111_ip_pushButton.clicked.connect(self.display_help_111_ip)
        self._mw.save_111_ip_pushButton.clicked.connect(self.save_routine)
        self._mw.Bx_111_ip_checkBox.stateChanged.connect(self.update_plot)
        self._mw.By_111_ip_checkBox.stateChanged.connect(self.update_plot)
        self._mw.Bz_111_ip_checkBox.stateChanged.connect(self.update_plot)
        self._mw.BNV_111_ip_checkBox.stateChanged.connect(self.update_plot)

        # disable save when no data computed yet
        self._mw.save_111_ip_pushButton.setEnabled(False)

        # set up plot
        self.xmin = self._mw.xmin_111_ip_doubleSpinBox.value()
        self.xmax = self._mw.xmax_111_ip_doubleSpinBox.value()
        self.nb_pts = self._mw.nb_pts_111_ip_spinBox.value()

        self.r_array = np.linspace(self.xmin, self.xmax, self.nb_pts)*1e-9
        self.data = pd.DataFrame({"Bx": np.zeros(self.nb_pts),
                                  "By": np.zeros(self.nb_pts),
                                  "Bz": np.zeros(self.nb_pts),
                                  "BNV": np.zeros(self.nb_pts)},
                                 index=self.r_array)
        self.Bx_plot_111_ip = pg.PlotDataItem(self.r_array, self.data["Bx"].values,
                                           pen=pg.mkPen(color=(239, 83, 80), width=2))
        self.By_plot_111_ip = pg.PlotDataItem(self.r_array, self.data["By"].values,
                                           pen=pg.mkPen(color=(70, 25, 136), width=2))         
        self.Bz_plot_111_ip = pg.PlotDataItem(self.r_array, self.data["Bz"].values,
                                           pen=pg.mkPen(color=(22, 173, 170), width=2))
        self.BNV_plot_111_ip = pg.PlotDataItem(self.r_array, self.data["BNV"].values,
                                            pen=pg.mkPen(color=(255, 165, 0), width=2))
                                             

        self._mw.type1_111_ip_graphicsView.addItem(self.Bx_plot_111_ip)
        self._mw.type1_111_ip_graphicsView.addItem(self.By_plot_111_ip)
        self._mw.type1_111_ip_graphicsView.addItem(self.Bz_plot_111_ip)
        self._mw.type1_111_ip_graphicsView.addItem(self.BNV_plot_111_ip)
        self._mw.type1_111_ip_graphicsView.setLabel("bottom", "r", units="m")
        self._mw.type1_111_ip_graphicsView.setLabel("left", "B", units="T")
        self._mw.type1_111_ip_graphicsView.setBackground("w")

        legend = pg.LegendItem(colCount=4)
        legend.addItem(self.Bx_plot_111_ip, name="  Bx")
        legend.addItem(self.By_plot_111_ip, name="  By")
        legend.addItem(self.Bz_plot_111_ip, name="  Bz")
        legend.addItem(self.BNV_plot_111_ip, name="  BNV")
        legend.setParentItem(self._mw.type1_111_ip_graphicsView.getPlotItem())
        legend.setPos(100, -10)
        
        return
    
    
    def set_up_and_calc_001t1(self):
        """
        Does the actual computation
        """
        self.cycloid_type = "001, type 1"
        
        self.xmin = self._mw.xmin_001t1_doubleSpinBox.value()
        self.xmax = self._mw.xmax_001t1_doubleSpinBox.value()
        self.nb_pts = self._mw.nb_pts_001t1_spinBox.value()

        self.r_array = np.linspace(self.xmin, self.xmax, self.nb_pts)
        self.data = pd.DataFrame({"Bx": np.zeros(self.nb_pts),
                                  "By": np.zeros(self.nb_pts),
                                  "Bz": np.zeros(self.nb_pts),
                                  "BNV": np.zeros(self.nb_pts)},
                                 index=self.r_array)
        self.r_array = self.r_array*1e-9

        period = self._mw.period_001t1_doubleSpinBox.value()*1e-9
        mDM = self._mw.mDM_001t1_doubleSpinBox.value()
        t = self._mw.t_001t1_doubleSpinBox.value()*1e-9
        dNV = self._mw.dNV_001t1_doubleSpinBox.value()*1e-9
        theta = self._mw.theta_001t1_doubleSpinBox.value()*np.pi/180
        phi = self._mw.phi_001t1_doubleSpinBox.value()*np.pi/180
        
        angle = self._mw.prof_001t1_doubleSpinBox.value()*np.pi/180
        x = self.r_array*np.cos(angle)
        y = self.r_array*np.sin(angle)

        self.params_dict={"Cycloid type": "Type 1",
                          "Crystal orientation": "(001), P along [111]",
                          "Propagation vector" : self._mw.k_001t1_comboBox.currentText(),
                          "Period (nm)": self._mw.period_001t1_doubleSpinBox.value(),
                          "mDM (µB)": self._mw.mDM_001t1_doubleSpinBox.value(),
                          "BFO thickness (nm)": self._mw.t_001t1_doubleSpinBox.value(),
                          "Tip  height dNV (nm)": self._mw.dNV_001t1_doubleSpinBox.value(),
                          "Tip polar angle θ (°)": self._mw.theta_001t1_doubleSpinBox.value(),
                          "Tip azimuthal angle φ (°)": self._mw.phi_001t1_doubleSpinBox.value(),
                          "Profile angle (°)": self._mw.prof_001t1_doubleSpinBox.value(),
                          "Number of points": self._mw.nb_pts_001t1_spinBox.value()}
        
        if self._mw.k_001t1_comboBox.currentText()=="k1 [1-10]":
            self.data["Bx"] = type1.Bx_k1(x, y, dNV, period, mDM, t/type1.a)
            self.data["By"] = type1.By_k1(x, y, dNV, period, mDM, t/type1.a)
            self.data["Bz"] = type1.Bz_k1(x, y, dNV, period, mDM, t/type1.a)

        elif self._mw.k_001t1_comboBox.currentText()=="k2 [-101]":
            self.data["Bx"] = type1.Bx_k2(x, y, dNV, period, mDM, t/type1.a)
            self.data["By"] = type1.By_k2(x, y, dNV, period, mDM, t/type1.a)
            self.data["Bz"] = type1.Bz_k2(x, y, dNV, period, mDM, t/type1.a)

        elif self._mw.k_001t1_comboBox.currentText()=="k3 [0-11]":
            self.data["Bx"] = type1.Bx_k3(x, y, dNV, period, mDM, t/type1.a)
            self.data["By"] = type1.By_k3(x, y, dNV, period, mDM, t/type1.a)
            self.data["Bz"] = type1.Bz_k3(x, y, dNV, period, mDM, t/type1.a)


        self.data["BNV"] = project_B(self.data["Bx"], self.data["By"], self.data["Bz"],
                                     theta, phi)

        self._mw.save_001t1_pushButton.setEnabled(True)
        self._mw.save_001t2_pushButton.setEnabled(False)
        self._mw.save_111_pushButton.setEnabled(False)
        self._mw.save_111_ip_pushButton.setEnabled(False)
        self.update_plot()
        return
    
    
    def set_up_and_calc_001t2(self):
        self.cycloid_type = "001, type 2"
        
        self.xmin = self._mw.xmin_001t2_doubleSpinBox.value()
        self.xmax = self._mw.xmax_001t2_doubleSpinBox.value()
        self.nb_pts = self._mw.nb_pts_001t2_spinBox.value()

        self.r_array = np.linspace(self.xmin, self.xmax, self.nb_pts)
        self.data = pd.DataFrame({"Bx": np.zeros(self.nb_pts),
                                  "By": np.zeros(self.nb_pts),
                                  "Bz": np.zeros(self.nb_pts),
                                  "BNV": np.zeros(self.nb_pts)},
                                 index=self.r_array)
        self.r_array = self.r_array*1e-9

        period = self._mw.period_001t2_doubleSpinBox.value()*1e-9
        mDM = self._mw.mDM_001t2_doubleSpinBox.value()
        t = self._mw.t_001t2_doubleSpinBox.value()*1e-9
        dNV = self._mw.dNV_001t2_doubleSpinBox.value()*1e-9
        theta = self._mw.theta_001t2_doubleSpinBox.value()*np.pi/180
        phi = self._mw.phi_001t2_doubleSpinBox.value()*np.pi/180
        
        angle = self._mw.prof_001t2_doubleSpinBox.value()*np.pi/180
        x = self.r_array*np.cos(angle)
        y = self.r_array*np.sin(angle)

        self.params_dict={"Cycloid type": "Type 2",
                          "Crystal orientation": "(001), P along [111]",
                          "Propagation vector" : self._mw.k_001t2_comboBox.currentText(),
                          "Period (nm)": self._mw.period_001t2_doubleSpinBox.value(),
                          "mDM (µB)": self._mw.mDM_001t2_doubleSpinBox.value(),
                          "BFO thickness (nm)": self._mw.t_001t2_doubleSpinBox.value(),
                          "Tip  height dNV (nm)": self._mw.dNV_001t2_doubleSpinBox.value(),
                          "Tip polar angle θ (°)": self._mw.theta_001t2_doubleSpinBox.value(),
                          "Tip azimuthal angle φ (°)": self._mw.phi_001t2_doubleSpinBox.value(),
                          "Profile angle (°)": self._mw.prof_001t2_doubleSpinBox.value(),
                          "Number of points": self._mw.nb_pts_001t2_spinBox.value()}
        
        if self._mw.k_001t2_comboBox.currentText()=="k'1 [-211]":
            self.data["Bx"] = type2.Bx_k1(x, y, dNV, period, mDM, t/type1.a)
            self.data["By"] = type2.By_k1(x, y, dNV, period, mDM, t/type1.a)
            self.data["Bz"] = type2.Bz_k1(x, y, dNV, period, mDM, t/type1.a)

        elif self._mw.k_001t2_comboBox.currentText()=="k'2 [1-21]":
            self.data["Bx"] = type2.Bx_k2(x, y, dNV, period, mDM, t/type1.a)
            self.data["By"] = type2.By_k2(x, y, dNV, period, mDM, t/type1.a)
            self.data["Bz"] = type2.Bz_k2(x, y, dNV, period, mDM, t/type1.a)

        elif self._mw.k_001t2_comboBox.currentText()=="k'3 [11-2]":
            alpha = self._mw.alpha_doubleSpinBox.value()*np.pi/180
            self.data["Bx"] = type2.Bx_k3_alpha(x, y, dNV, period, mDM, t/type1.a, alpha)
            self.data["By"] = type2.By_k3_alpha(x, y, dNV, period, mDM, t/type1.a, alpha)
            self.data["Bz"] = type2.Bz_k3_alpha(x, y, dNV, period, mDM, t/type1.a, alpha)


        self.data["BNV"] = project_B(self.data["Bx"], self.data["By"], self.data["Bz"],
                                     theta, phi)

        self._mw.save_001t1_pushButton.setEnabled(False)
        self._mw.save_001t2_pushButton.setEnabled(True)
        self._mw.save_111_pushButton.setEnabled(False)
        self._mw.save_111_ip_pushButton.setEnabled(False)
        self.update_plot()
        return

    
    def set_up_and_calc_111(self):
        self.cycloid_type = "111"

        self.xmin = self._mw.xmin_111_doubleSpinBox.value()
        self.xmax = self._mw.xmax_111_doubleSpinBox.value()
        self.nb_pts = self._mw.nb_pts_111_spinBox.value()

        self.r_array = np.linspace(self.xmin, self.xmax, self.nb_pts)
        self.data = pd.DataFrame({"Bx": np.zeros(self.nb_pts),
                                  "By": np.zeros(self.nb_pts),
                                  "Bz": np.zeros(self.nb_pts),
                                  "BNV": np.zeros(self.nb_pts)},
                                 index=self.r_array)
        self.r_array = self.r_array*1e-9

        period = self._mw.period_111_doubleSpinBox.value()*1e-9
        mS = self._mw.mS_111_doubleSpinBox.value()
        N = self._mw.t_111_spinBox.value()
        dNV = self._mw.dNV_111_doubleSpinBox.value()*1e-9
        theta = self._mw.theta_111_doubleSpinBox.value()*np.pi/180
        phi = self._mw.phi_111_doubleSpinBox.value()*np.pi/180

        angle = self._mw.prof_111_doubleSpinBox.value()*np.pi/180
        x = self.r_array*np.cos(angle)
        y = self.r_array*np.sin(angle)

        self.params_dict={"Cycloid type": "Type 1",
                          "Crystal orientation": "(111), P along [111]",
                          "Period (nm)": self._mw.period_111_doubleSpinBox.value(),
                          "mS (µB)": self._mw.mS_111_doubleSpinBox.value(),
                          "Number of BFO layers": self._mw.t_111_spinBox.value(),
                          "Tip  height dNV (nm)": self._mw.dNV_111_doubleSpinBox.value(),
                          "Tip polar angle θ (°)": self._mw.theta_111_doubleSpinBox.value(),
                          "Tip azimuthal angle φ (°)": self._mw.phi_111_doubleSpinBox.value(),
                          "Profile angle (°)": self._mw.prof_111_doubleSpinBox.value(),
                          "Number of points": self._mw.nb_pts_111_spinBox.value()}

        self.data["Bx"] = t1_111.Bx(x, y, dNV, period, mS, N)
        self.data["By"] = t1_111.By(x, y, dNV, period, mS, N)
        self.data["Bz"] = t1_111.Bz(x, y, dNV, period, mS, N)

        self.data["BNV"] = project_B(self.data["Bx"], self.data["By"], self.data["Bz"],
                                     theta, phi)

        self._mw.save_001t1_pushButton.setEnabled(False)
        self._mw.save_001t2_pushButton.setEnabled(False)
        self._mw.save_111_pushButton.setEnabled(True)
        self._mw.save_111_ip_pushButton.setEnabled(False)
        self.update_plot()
        return


    def set_up_and_calc_111_ip(self): 
        self.cycloid_type = "IP 111"

        self.xmin = self._mw.xmin_111_ip_doubleSpinBox.value()
        self.xmax = self._mw.xmax_111_ip_doubleSpinBox.value()
        self.nb_pts = self._mw.nb_pts_111_ip_spinBox.value()

        self.r_array = np.linspace(self.xmin, self.xmax, self.nb_pts)
        self.data = pd.DataFrame({"Bx": np.zeros(self.nb_pts),
                                  "By": np.zeros(self.nb_pts),
                                  "Bz": np.zeros(self.nb_pts),
                                  "BNV": np.zeros(self.nb_pts)},
                                 index=self.r_array)
        self.r_array = self.r_array*1e-9

        period = self._mw.period_111_ip_doubleSpinBox.value()*1e-9
        mS = self._mw.mS_111_ip_doubleSpinBox.value()
        mDM = self._mw.mDM_111_ip_doubleSpinBox.value()
        phase = self._mw.phase_111_ip_doubleSpinBox.value()
        N = self._mw.t_111_ip_spinBox.value()
        dNV = self._mw.dNV_111_ip_doubleSpinBox.value()*1e-9
        theta = self._mw.theta_111_ip_doubleSpinBox.value()*np.pi/180
        phi = self._mw.phi_111_ip_doubleSpinBox.value()*np.pi/180
        if self._mw.rot_sense_111_ip_comboBox.currentText() == "CW":
            eps = -1
        else:
            eps = 1

        angle = self._mw.prof_111_ip_doubleSpinBox.value()*np.pi/180
        x = self.r_array*np.cos(angle)
        y = self.r_array*np.sin(angle)

        self.params_dict={"Cycloid type": "Type 1",
                          "Crystal orientation": "(111), P along [11-1]",
                          "Propagation vector" : self._mw.k_111_ip_comboBox.currentText(),
                          "Period (nm)": self._mw.period_111_ip_doubleSpinBox.value(),
                          "mS (µB)": self._mw.mS_111_ip_doubleSpinBox.value(),
                          "mDM (µB)": self._mw.mDM_111_ip_doubleSpinBox.value(),
                          "Phase (°)": self._mw.phase_111_ip_doubleSpinBox.value(),
                          "Rotational sense": self._mw.rot_sense_111_ip_comboBox.currentText(),
                          "Number of BFO layers": self._mw.t_111_ip_spinBox.value(),
                          "Tip  height dNV (nm)": self._mw.dNV_111_ip_doubleSpinBox.value(),
                          "Tip polar angle θ (°)": self._mw.theta_111_ip_doubleSpinBox.value(),
                          "Tip azimuthal angle φ (°)": self._mw.phi_111_ip_doubleSpinBox.value(),
                          "Profile angle (°)": self._mw.prof_111_ip_doubleSpinBox.value(),
                          "Number of points": self._mw.nb_pts_111_ip_spinBox.value()}
        
        if self._mw.k_111_ip_comboBox.currentText()=="k1":
            self.data["Bx"] = t1_111_ip.Bx_k1(x, y, dNV, period, mS, mDM, N, eps, phase)
            self.data["By"] = t1_111_ip.By_k1(x, y, dNV, period, mS, mDM, N, eps, phase)
            self.data["Bz"] = t1_111_ip.Bz_k1(x, y, dNV, period, mS, mDM, N, eps, phase)

        elif self._mw.k_111_ip_comboBox.currentText()=="k2":
            self.data["Bx"] = t1_111_ip.Bx_k2(x, y, dNV, period, mS, mDM, N, eps, phase)
            self.data["By"] = t1_111_ip.By_k2(x, y, dNV, period, mS, mDM, N, eps, phase)
            self.data["Bz"] = t1_111_ip.Bz_k2(x, y, dNV, period, mS, mDM, N, eps, phase)

        elif self._mw.k_111_ip_comboBox.currentText()=="k3":
            self.data["Bx"] = t1_111_ip.Bx_k3(x, y, dNV, period, mS, mDM, N, eps, phase)
            self.data["By"] = t1_111_ip.By_k3(x, y, dNV, period, mS, mDM, N, eps, phase)
            self.data["Bz"] = t1_111_ip.Bz_k3(x, y, dNV, period, mS, mDM, N, eps, phase)

        self.data["BNV"] = project_B(self.data["Bx"], self.data["By"], self.data["Bz"],
                                     theta, phi)

        self._mw.save_001t1_pushButton.setEnabled(False)
        self._mw.save_001t2_pushButton.setEnabled(False)
        self._mw.save_111_pushButton.setEnabled(False)
        self._mw.save_111_ip_pushButton.setEnabled(True)
        self.update_plot()
        return
    

    def update_plot(self):
        """ Plots the data. """
        
        if np.abs(np.mean(np.gradient(self.data["BNV"]))) > 1e-10:
            minima = argrelmin(np.abs(self.data["BNV"].values))[0]
            if len(minima) >=2:
                est_period = 2*(self.data.index.values[minima[1]]\
                                -self.data.index.values[minima[0]])
            else:
                est_period = 2*(self.data.index[-1]-self.data.index[0])
            p0=[np.abs(np.max(self.data["BNV"])), est_period, 0.5]
            popt, pcov = curve_fit(sinus, self.data.index, self.data["BNV"].values, p0=p0)
            periodic = True
        else:
            periodic = False
        
        if self.cycloid_type == "001, type 1":
            self.Bx_plot_001t1.setData(x=self.r_array, y=self.data["Bx"])
            self.By_plot_001t1.setData(x=self.r_array, y=self.data["By"])
            self.Bz_plot_001t1.setData(x=self.r_array, y=self.data["Bz"])
            self.BNV_plot_001t1.setData(x=self.r_array, y=self.data["BNV"])

            plot_list = self._mw.type1_001_graphicsView.listDataItems()
            
            if self._mw.Bx_001t1_checkBox.isChecked():
                if self.Bx_plot_001t1 not in plot_list:
                    self._mw.type1_001_graphicsView.addItem(self.Bx_plot_001t1)
            else:
                if self.Bx_plot_001t1 in plot_list:
                    self._mw.type1_001_graphicsView.removeItem(self.Bx_plot_001t1)
                
            if self._mw.By_001t1_checkBox.isChecked():
                if self.By_plot_001t1 not in plot_list:
                    self._mw.type1_001_graphicsView.addItem(self.By_plot_001t1)
            else:
                if self.By_plot_001t1 in plot_list:
                    self._mw.type1_001_graphicsView.removeItem(self.By_plot_001t1)
                
            if self._mw.Bz_001t1_checkBox.isChecked():
                if self.Bz_plot_001t1 not in plot_list:
                    self._mw.type1_001_graphicsView.addItem(self.Bz_plot_001t1)
            else:
                if self.Bz_plot_001t1 in plot_list:
                    self._mw.type1_001_graphicsView.removeItem(self.Bz_plot_001t1)

            if self._mw.BNV_001t1_checkBox.isChecked():
                if self.BNV_plot_001t1 not in plot_list:
                    self._mw.type1_001_graphicsView.addItem(self.BNV_plot_001t1)
            else:
                if self.BNV_plot_001t1 in plot_list:
                    self._mw.type1_001_graphicsView.removeItem(self.BNV_plot_001t1)

            if periodic:
                self._mw.period_001t1_label.setText("Apparent period: {:.1f} nm".format(popt[1]))
                self._mw.ampl_001t1_label.setText("B_NV amplitude: {:.1f} µT".format(
                    np.abs(popt[0])*1e6))
            else:
                self._mw.period_001t1_label.setText("Apparent period: - ")
                self._mw.ampl_001t1_label.setText("B_NV amplitude: {:.1f} µT".format(
                    np.mean(self.data["BNV"]*1e6)))

        elif self.cycloid_type == "001, type 2":
            self.Bx_plot_001t2.setData(x=self.r_array, y=self.data["Bx"])
            self.By_plot_001t2.setData(x=self.r_array, y=self.data["By"])
            self.Bz_plot_001t2.setData(x=self.r_array, y=self.data["Bz"])
            self.BNV_plot_001t2.setData(x=self.r_array, y=self.data["BNV"])

            plot_list = self._mw.type2_001_graphicsView.listDataItems()
            
            if self._mw.Bx_001t2_checkBox.isChecked():
                if self.Bx_plot_001t2 not in plot_list:
                    self._mw.type2_001_graphicsView.addItem(self.Bx_plot_001t2)
            else:
                if self.Bx_plot_001t2 in plot_list:
                    self._mw.type2_001_graphicsView.removeItem(self.Bx_plot_001t2)
                
            if self._mw.By_001t2_checkBox.isChecked():
                if self.By_plot_001t2 not in plot_list:
                    self._mw.type2_001_graphicsView.addItem(self.By_plot_001t2)
            else:
                if self.By_plot_001t2 in plot_list:
                    self._mw.type2_001_graphicsView.removeItem(self.By_plot_001t2)
                
            if self._mw.Bz_001t2_checkBox.isChecked():
                if self.Bz_plot_001t2 not in plot_list:
                    self._mw.type2_001_graphicsView.addItem(self.Bz_plot_001t2)
            else:
                if self.Bz_plot_001t2 in plot_list:
                    self._mw.type2_001_graphicsView.removeItem(self.Bz_plot_001t2)

            if self._mw.BNV_001t2_checkBox.isChecked():
                if self.BNV_plot_001t2 not in plot_list:
                    self._mw.type2_001_graphicsView.addItem(self.BNV_plot_001t2)
            else:
                if self.BNV_plot_001t2 in plot_list:
                    self._mw.type2_001_graphicsView.removeItem(self.BNV_plot_001t2)

            if periodic:
                self._mw.period_001t2_label.setText("Apparent period: {:.1f} nm".format(popt[1]))
                self._mw.ampl_001t2_label.setText("B_NV amplitude: {:.1f} µT".format(
                    np.abs(popt[0])*1e6))
            else:
                self._mw.period_001t2_label.setText("Apparent period: - ")
                self._mw.ampl_001t2_label.setText("B_NV amplitude: {:.1f} µT".format(
                    np.mean(self.data["BNV"]*1e6)))

        elif self.cycloid_type == "111":
            self.Bx_plot_111.setData(x=self.r_array, y=self.data["Bx"])
            self.By_plot_111.setData(x=self.r_array, y=self.data["By"])
            self.Bz_plot_111.setData(x=self.r_array, y=self.data["Bz"])
            self.BNV_plot_111.setData(x=self.r_array, y=self.data["BNV"])

            plot_list = self._mw.type1_111_graphicsView.listDataItems()
            
            if self._mw.Bx_111_checkBox.isChecked():
                if self.Bx_plot_111 not in plot_list:
                    self._mw.type1_111_graphicsView.addItem(self.Bx_plot_111)
            else:
                if self.Bx_plot_111 in plot_list:
                    self._mw.type1_111_graphicsView.removeItem(self.Bx_plot_111)
                
            if self._mw.By_111_checkBox.isChecked():
                if self.By_plot_111 not in plot_list:
                    self._mw.type1_111_graphicsView.addItem(self.By_plot_111)
            else:
                if self.By_plot_111 in plot_list:
                    self._mw.type1_111_graphicsView.removeItem(self.By_plot_111)
                
            if self._mw.Bz_111_checkBox.isChecked():
                if self.Bz_plot_111 not in plot_list:
                    self._mw.type1_111_graphicsView.addItem(self.Bz_plot_111)
            else:
                if self.Bz_plot_111 in plot_list:
                    self._mw.type1_111_graphicsView.removeItem(self.Bz_plot_111)

            if self._mw.BNV_111_checkBox.isChecked():
                if self.BNV_plot_111 not in plot_list:
                    self._mw.type1_111_graphicsView.addItem(self.BNV_plot_111)
            else:
                if self.BNV_plot_111 in plot_list:
                    self._mw.type1_111_graphicsView.removeItem(self.BNV_plot_111)

            if periodic:
                self._mw.period_111_label.setText("Apparent period: {:.1f} nm".format(popt[1]))
                self._mw.ampl_111_label.setText("B_NV amplitude: {:.1f} µT".format(
                    np.abs(popt[0])*1e6))
            else:
                self._mw.period_111_label.setText("Apparent period: - ")
                self._mw.ampl_111_label.setText("B_NV amplitude: {:.1f} µT".format(
                    np.mean(self.data["BNV"]*1e6)))

        elif self.cycloid_type == "IP 111":
            self.Bx_plot_111_ip.setData(x=self.r_array, y=self.data["Bx"])
            self.By_plot_111_ip.setData(x=self.r_array, y=self.data["By"])
            self.Bz_plot_111_ip.setData(x=self.r_array, y=self.data["Bz"])
            self.BNV_plot_111_ip.setData(x=self.r_array, y=self.data["BNV"])

            plot_list = self._mw.type1_111_ip_graphicsView.listDataItems()
            
            if self._mw.Bx_111_ip_checkBox.isChecked():
                if self.Bx_plot_111_ip not in plot_list:
                    self._mw.type1_111_ip_graphicsView.addItem(self.Bx_plot_111_ip)
            else:
                if self.Bx_plot_111_ip in plot_list:
                    self._mw.type1_111_ip_graphicsView.removeItem(self.Bx_plot_111_ip)
                
            if self._mw.By_111_ip_checkBox.isChecked():
                if self.By_plot_111_ip not in plot_list:
                    self._mw.type1_111_ip_graphicsView.addItem(self.By_plot_111_ip)
            else:
                if self.By_plot_111_ip in plot_list:
                    self._mw.type1_111_ip_graphicsView.removeItem(self.By_plot_111_ip)
                
            if self._mw.Bz_111_ip_checkBox.isChecked():
                if self.Bz_plot_111_ip not in plot_list:
                    self._mw.type1_111_ip_graphicsView.addItem(self.Bz_plot_111_ip)
            else:
                if self.Bz_plot_111_ip in plot_list:
                    self._mw.type1_111_ip_graphicsView.removeItem(self.Bz_plot_111_ip)

            if self._mw.BNV_111_ip_checkBox.isChecked():
                if self.BNV_plot_111_ip not in plot_list:
                    self._mw.type1_111_ip_graphicsView.addItem(self.BNV_plot_111_ip)
            else:
                if self.BNV_plot_111_ip in plot_list:
                    self._mw.type1_111_ip_graphicsView.removeItem(self.BNV_plot_111_ip)

            if periodic:
                self._mw.period_111_ip_label.setText("Apparent period: {:.1f} nm".format(popt[1]))
                self._mw.ampl_111_ip_label.setText("B_NV amplitude: {:.1f} µT".format(
                    np.abs(popt[0])*1e6))
            else:
                self._mw.period_111_ip_label.setText("Apparent period: - ")
                self._mw.ampl_111_ip_label.setText("B_NV amplitude: {:.1f} µT".format(
                    np.mean(self.data["BNV"]*1e6)))
            
        return
                                         
        
    def save_routine(self):
        self.filedialog = QtWidgets.QFileDialog()
        self.filedialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
        self.textfile = self.filedialog.getSaveFileName(self._mw, "Choose a file name", "", "Text Files (*.txt)")
        self.data_to_save = pd.DataFrame({"Bx(µT)": self.data["Bx"]*1e6,
                                          "By(µT)": self.data["By"]*1e6,
                                          "Bz(µT)": self.data["Bz"]*1e6,
                                          "B_NV(µT)": self.data["BNV"]*1e6},
                                         index=self.data.index)
        header = "======\nHeader\n======\n"+pp.pformat(self.params_dict,
                                                       indent=0)[1:-1]+"\n\n======\nData\n======\n"
        with open(self.textfile[0], "w") as f:
            f.write(header)
        self.data_to_save.to_csv(self.textfile[0], sep=" ", index_label="r(nm)",
                                 quotechar=" ", header=True, mode="a")


    def display_help_001t1(self):
        self.h_001t1_dlg = QtWidgets.QDialog()
        self.h_001t1_dlg.setWindowTitle("Help for cycloid type 1, BFO 001")
        text_label = QtWidgets.QLabel(" ", self.h_001t1_dlg)
        text_label.setText("In this configuration, the surface of the BFO film is parallel to\n"
                           "the (001) plane and the ferroelectric polarization direction is [111].\n \n"
                           "We are interested in the cycloid type 1, with the\n"
                           "propagation vectors along three possible directions:\n"
                           "  - k1, parallel to [1-10]\n"
                           "  - k2, parallel to [-101]\n"
                           "  - k3, parallel to [0-11]\n\n"
                           "The parameter 'Period' is the bulk period of the cycloid,\n"
                           "usually 64 nm.\n"
                           "m_DM is the strength of the DMI induced magnetic moments, in µB.\n"
                           "These moments are assumed to be in the plane perpendicular to P\n" 
                           "and k.\n"
                           "The tip angle θ is the polar angle, θ = 0° means along [001], \n"
                           "and φ is the azimuthal angle, measured from the [100] direction.\n"
                           "The profile angle is also measured with respect to [100].\n"
                           )
        text_label.move(20, 20)
        button = QtWidgets.QPushButton("OK", self.h_001t1_dlg)
        button.clicked.connect(self.h_001t1_dlg.accept)
        button.move(305, 320)
        self.h_001t1_dlg.exec()


    def display_help_001t2(self):
        self.h_001t2_dlg = QtWidgets.QDialog()
        self.h_001t2_dlg.setWindowTitle("Help for cycloid type 2, BFO 001")
        text_label = QtWidgets.QLabel(" ", self.h_001t2_dlg)
        text_label.setText("In this configuration, the surface of the BFO film is parallel to\n"
                           "the (001) plane and the ferroelectric polarization direction is [111].\n \n"
                           "We are interested in the cycloid type 2, with the\n"
                           "propagation vectors along three possible directions:\n"
                           "  - k'1, parallel to [-211]\n"
                           "  - k'2, parallel to [1-21]\n"
                           "  - k'3, parallel to [11-2]\n\n"
                           "In the case of k'3, no stray field is expected in the ideal configuration.\n"
                           "You need therefore to add a deviation angle α, corresponding to a rotation\n"
                           "of k'3 in the plane perpendicular to P.\n \n"
                           "The parameter 'Period' is the bulk period of the cycloid.\n"
                           "m_DM is the strength of the DMI induced magnetic moments, in µB.\n"
                           "These moments are assumed to be in the plane perpendicular to P and k.\n"
                           "The tip angle θ is the polar angle, θ=0° means along [001], \n"
                           "and φ is the azimuthal angle, measured from the [100] direction.\n"
                           "The profile angle is also measured with respect to [100].\n"
                           )
        text_label.move(20, 20)
        button = QtWidgets.QPushButton("OK", self.h_001t2_dlg)
        button.clicked.connect(self.h_001t2_dlg.accept)
        button.move(350, 355)
        self.h_001t2_dlg.exec()


    def display_help_111(self):
        self.h_111_dlg = QtWidgets.QDialog()
        self.h_111_dlg.setWindowTitle("Help for cycloid type 1, BFO 111, P out-of-plane")
        text_label = QtWidgets.QLabel(" ", self.h_111_dlg)
        text_label.setText("In this configuration, the surface of the BFO film is parallel to\n"
                           "the (111) plane and the ferroelectric polarization direction is [111],\n"
                           "meaning out-of-plane\n \n"
                           "We are interested in the cycloid type 1, and the three propagation\n"
                           "directions are in-plane and equivalent.\n"
                           "k is chosen along x, so By is always 0.\n \n"
                           "The parameter 'Period' is the bulk period of the cycloid, usually 64 nm.\n"
                           "m_S is the strength of the Fe magnetic moments, in µB, and\n"
                           "the DMI induced moments do not contribute in this configuration.\n \n"
                           "The number of BFO layers is critical in this case, an odd number of\n"
                           "layers produces a much stronger field than an even number because\n"
                           "the last layer is not compensated.\n \n"
                           "The tip angle θ is the polar angle, θ=0° means along [111], \n"
                           "and φ is the azimuthal angle, measured from the k direction.\n"
                           "The profile angle is also measured with respect to k.\n"
                           )
        text_label.move(20, 20)
        button = QtWidgets.QPushButton("OK", self.h_111_dlg)
        button.clicked.connect(self.h_111_dlg.accept)
        button.move(330, 340)
        self.h_111_dlg.exec()


    def display_help_111_ip(self):
        self.h_111_ip_dlg = QtWidgets.QDialog()
        self.h_111_ip_dlg.setWindowTitle("Help for cycloid type 1, BFO 111, P mostly in-plane")
        text_label = QtWidgets.QLabel(" ", self.h_111_ip_dlg)
        text_label.setText("In this configuration, the surface of the BFO film is parallel to\n"
                           "the (111) plane and the ferroelectric polarization direction is [11-1],\n"
                           "meaning out-of-plane\n \n"
                           "We are interested in the cycloid type 1, with the\n"
                           "propagation vectors along three possible directions:\n"
                           "  - k1, parallel to [1-10], corresponding to x in the lab frame\n"
                           "  - k2, parallel to [011]\n"
                           "  - k3, parallel to [101]\n\n"
                           "The parameter 'Period' is the bulk period of the cycloid, usually 64 nm.\n"
                           "m_S is the strength of the Fe magnetic moments, in µB.\n"
                           "m_DM is the strength of the DMI induced magnetic moments, in µB.\n"
                           "The 'Phase' parameter corresponds to the dephasing (in degrees)\n"
                           "between the cycloid and the spin density wave. \n"
                           "The rotational sense of the cycloid plays also a role, and going\n"
                           "from P to -P amounts to a change of sense.\n"
                           "The number of BFO layers is critical in this case, an odd number of\n"
                           "layers produces a stronger field than an even number because\n"
                           "the last layer is not compensated.\n \n"
                           "The tip angle θ is the polar angle, θ=0° means along [111], \n"
                           "and φ is the azimuthal angle, measured from the k1 direction.\n"
                           "The profile angle is also measured with respect to k.\n"
                           )
        text_label.move(20, 20)
        button = QtWidgets.QPushButton("OK", self.h_111_ip_dlg)
        button.clicked.connect(self.h_111_ip_dlg.accept)
        button.move(330, 450)
        self.h_111_ip_dlg.exec()
        

def main():
    # Create an application object.
    app = QtWidgets.QApplication(sys.argv)
    gui = BFOator()
    gui.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
