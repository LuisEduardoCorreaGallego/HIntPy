#! /usr/bin/env python
# -*- coding: utf-8 -*-
u"""
Name:
File Name: PinchAnalysis
Description: None
Category: Chemical Engineering, Heat Integration
Requested Elements:
Author: Luis Eduardo Correa Gallego <luise.correa@udea.edu.co>
Created on: 5/01/2019
Last modification: 6/01/2019
Used IDE: PyCharm Professional Edition
"""
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import numpy as np
import os
from PIL import Image, ImageTk
import Pmw
from tkinter import ttk
import tkinter
from tkinter.scrolledtext import ScrolledText
matplotlib.use('TkAgg')
plt.style.use('seaborn')


class pinchStream(object):
    """
    This class
    """
    def __init__(self, file):
        self.file = str(file)
        global Specifications_array, Specifications
        global T_in, T_out, Q_dot
        global heatFlowSink, ptCascade, temperatures
        global tempSink, tempSource, heatFlowSource, heatFlowSink, CP_stream

        # Load data from text file
        Specifications_array = np.loadtxt(self.file,
                                          skiprows=1,
                                          unpack=False)
        # Generate stream list from information on text file
        streamList = {}
        for i in range(0, len(Specifications_array)):
            streamList["Stream_{0}".format(i+1)] = list(Specifications_array[i])
        locals().update(streamList)
        Specifications = list(streamList.values())
        # Separating data by type: heat flow,
        #                          start temperature,
        #                          target temperature
        #                          approach
        Data = [[Specifications[i][0] for i in range(len(Specifications))],
                [Specifications[i][1] for i in range(len(Specifications))],
                [Specifications[i][2] for i in range(len(Specifications))],
                [Specifications[i][3] for i in range(len(Specifications))]]
        # Generate sets for each type
        Q_dot = Data[0]
        T_in = Data[1]
        T_out = Data[2]
        dT_min = Data[3]
        # Calculate the CP for each stream
        CP_stream = np.asarray(Q_dot) / (np.asarray(T_out)-np.asarray(T_in))
        # Get hot and cold stream index
        hotStreamIndex = []
        coldStreamIndex = []
        for i in range(len(Specifications)):
            if Specifications[i][1] > Specifications[i][2]:
                hotStreamIndex.append(i)
            elif Specifications[i][1] < Specifications[i][2]:
                coldStreamIndex.append(i)
        # Calculate source and sink heat flow
        Q_dot_source = 0
        for i in range(len(hotStreamIndex)):
            Q_dot_source = Q_dot_source + Specifications[hotStreamIndex[i]][0]
        Q_dot_sink = 0
        for i in range(len(coldStreamIndex)):
            Q_dot_sink = Q_dot_sink + Specifications[coldStreamIndex[i]][0]
        # Shifted temperature and heat capacity flow hot stream
        for i in range(len(hotStreamIndex)):
            deltaT = 0.5 * Specifications[i][3]
            row = hotStreamIndex[i]
            Specifications[row].append(Specifications[row][1] - deltaT)
            Specifications[row].append(Specifications[row][2] - deltaT)
            Specifications[row].append(Specifications[row][0] / (
                        Specifications[row][5] - Specifications[row][4]))
        # Shifted temperature and heat capacity flow cold stream
        deltaT = 0
        for i in range(len(coldStreamIndex)):
            deltaT = 0.5 * Specifications[i][3]
            row = coldStreamIndex[i]
            Specifications[row].append(Specifications[row][1] + deltaT)
            Specifications[row].append(Specifications[row][2] + deltaT)
            Specifications[row].append(Specifications[row][0] / (
                        Specifications[row][5] - Specifications[row][4]))
        # Get temperatures and intervals
        temperatures = []
        for i in range(len(Specifications)):
            temperatures.append(Specifications[i][4])
            temperatures.append(Specifications[i][5])
        # Get sorting index
        tempInd = np.argsort(temperatures)
        # Get rid of duplicates
        temperatures = set(temperatures)
        temperatures = list(temperatures)
        # Sort temperatures
        temperatures.sort(reverse=True)
        # Set cascade temperatures
        cascade = []
        cascadeSink = []
        cascadeSource = []
        for i in range(len(temperatures) - 1):
            cascade.append([])
            cascadeSink.append([])
            cascadeSource.append([])
            cascade[i].append(temperatures[i])
            cascadeSink[i].append(temperatures[i])
            cascadeSource[i].append(temperatures[i])
            cascade[i].append(temperatures[i + 1])
            cascadeSink[i].append(temperatures[i + 1])
            cascadeSource[i].append(temperatures[i + 1])
        # ===
        for i in range(len(cascade)):
            cascade[i].append(cascade[i][0] - cascade[i][1])
            cascadeSink[i].append(cascade[i][0] - cascade[i][1])
            cascadeSource[i].append(cascade[i][0] - cascade[i][1])
            C_source = 0
            C_sink = 0
            for j in range(len(Specifications)):
                if (cascade[i][0] <= Specifications[j][4]) and (
                        cascade[i][1] >= Specifications[j][5]) and (
                        Specifications[j][6] < 0):
                    C_source = C_source + Specifications[j][6]
            for j in range(len(Specifications)):
                if (cascade[i][0] > Specifications[j][4]) and (
                        cascade[i][1] < Specifications[j][5]) and (
                        Specifications[j][6] > 0):
                    C_sink = C_sink + Specifications[j][6]
            cascade[i].append(C_source + C_sink)
            cascade[i].append(cascade[i][3] * cascade[i][2])
            cascadeSink[i].append(C_sink)
            cascadeSink[i].append(cascadeSink[i][3] * cascadeSink[i][2])
            cascadeSource[i].append(C_source)
            cascadeSource[i].append(cascadeSource[i][3] * cascadeSource[i][2])
        heatFlow = 0
        ptCascade = []
        ptCascade.append(heatFlow)
        for i in range(len(cascade)):
            heatFlow = heatFlow - cascade[i][4]
            ptCascade.append(heatFlow)
        # Calculate recovered and remaining heat flow
        if [n for n in ptCascade if n < 0]:
            Q_dot_rem = min([n for n in ptCascade if n < 0])
        else:
            Q_dot_rem = 0
        Q_dot_rec = Q_dot_sink - abs(Q_dot_rem)
        heatFlow = abs(Q_dot_rem)
        ptCascade = []
        ptCascade.append(heatFlow)
        for i in range(len(cascade)):
            heatFlow = heatFlow - cascade[i][4]
            ptCascade.append(heatFlow)
        # Prepare plot for source
        heatFlowSource = []
        tempSource = []
        heatFlowSource.append(0.)
        for i in reversed(cascadeSource):
            if i[4] < 0:
                heatFlowSource.append(heatFlowSource[-1] + i[4] * -1.)
                tempSource.append(i[0] + deltaT)
                tempSource.append(i[1] + deltaT)
        # Get rid of duplicates
        tempSource = set(tempSource)
        tempSource = list(tempSource)
        # Sort temperatures
        tempSource.sort(reverse=False)
        # Prepare plot for sink
        heatFlowSink = []
        tempSink = []
        # Generating heat flow cascade
        # The last element of cascade corresponds to start of sink
        heatFlowSink.append(ptCascade[-1])
        for i in reversed(cascadeSink):
            if i[4] > 0:
                heatFlowSink.append(heatFlowSink[-1] + i[4])
                tempSink.append(i[0] - deltaT)
                tempSink.append(i[1] - deltaT)
        # Get rid of duplicates
        tempSink = set(tempSink)
        tempSink = list(tempSink)
        # Sort temperatures
        tempSink.sort(reverse=False)

    @staticmethod
    def streamData(T_units='°C', Q_units='kW'):
        """
        This method shows the stream data with its information
        :param T_units: Units for temperature, default is '°C'
        :param Q_units: Units for heat, default is 'kW'
        :return: Table with stream specification, generated on tkinter
        """
        header = ("Stream", "Type",
                  "Inlet temperature \n ("+T_units+")",
                  "Outlet temperature \n ("+T_units+")",
                  "Available heat \n ("+Q_units+")")
        colors = ("brown4", "midnight blue", "midnight blue",
                  "midnight blue", "midnight blue")
        for Header, Color, Cell in zip(header, colors, range(len(header))):
            tkinter.Label(streamsSpecificationsTab, text=Header,
                          background=Color,
                          foreground="white",
                          relief="raised",
                          borderwidth=3).grid(row=0, column=Cell + 1,
                                              columnspan=1,
                                              rowspan=1,
                                              sticky='nsew')
        for i in range(len(T_in)):
            tkinter.Label(streamsSpecificationsTab, text=i + 1,
                          background="slate gray",
                          foreground="black",
                          relief="raised",
                          borderwidth=3).grid(row=i + 1, column=1,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for m in range(len(T_in)):
            if T_in[m] > T_out[m]:
                text = "Hot"
                color = "red"
            else:
                text = "Cold"
                color = "blue"
            tkinter.Label(streamsSpecificationsTab, text=text,
                          background="khaki1",
                          foreground=color,
                          relief="ridge",
                          borderwidth=3).grid(row=m + 1, column=2,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for j in range(len(T_in)):
            tkinter.Label(streamsSpecificationsTab, text=T_in[j],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=j + 1, column=3,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for k in range(len(T_in)):
            tkinter.Label(streamsSpecificationsTab, text=T_out[k],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=k + 1, column=4,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for n in range(len(T_in)):
            tkinter.Label(streamsSpecificationsTab, text=Q_dot[n],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=n + 1, column=5,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')

    @staticmethod
    def drawStreams(Data='1', T_units='°C', Q_units='kW'):
        """
        This method draws a plot of corresponding streams
        :param Data: Name for the set of streams
        :param T_units: Units for temperature, default is '°C'
        :param Q_units: Units for heat, default is 'kW'
        :return:
        """
        global figStreams, axStreams
        figStreams = plt.figure()
        figStreams.canvas.set_window_title('Streams for data')
        axStreams = figStreams.add_subplot(111)
        axStreams.minorticks_on()
        axStreams.grid(True, which='BOTH', ls=':')
        labels = []
        if Q_units == 'kW':
            Q_denom = 'kW'
            Q_numer = ''
        elif Q_units == 'BTU/h':
            Q_denom = 'BTU'
            Q_numer = 'h \\times'
        else:
            pass
        for j in range(len(T_in)+1):
            labels.append("Stream " +
                          str(j)+"\n" +
                          "$CP_{" + str(j)+"}" + "\ = \ $" +
                          str(round(abs(CP_stream[j-1]), 4)) +
                          "$ \ \\frac{" + Q_denom +
                          "}{" +
                          Q_numer +
                          T_units +
                          "}$")
        axStreams.set_yticklabels(labels)
        axStreams.set_xticklabels([])
        for k in range(len(T_in)):
            if T_in[k] > T_out[k]:
                colorLine = 'red'
            else:
                colorLine = 'blue'
            axStreams.plot([T_in[k], T_out[k]], [k+1, k+1], color=colorLine,
                           linestyle='-', linewidth=3, marker='o')
            axStreams.arrow(T_in[k], k+1, (T_out[k]-T_in[k])/2, 0,
                            head_width=0.2, head_length=4, color=colorLine)
            axStreams.annotate(str(T_in[k])+" "+T_units,
                               xy=(T_in[k], k+0.92), ha='center')
            axStreams.annotate(str(T_out[k])+" "+T_units,
                               xy=(T_out[k], k+0.92), ha='center')
        axStreams.yaxis.set_major_locator(MaxNLocator(integer=True))
        axStreams.invert_yaxis()
        figStreams.tight_layout()
        return figStreams.savefig('streamsCurve.jpg', bbox_inches="tight",
                                  metadata={'Author': 'Luis Eduardo Correa Gallego'},
                                  dpi=400)

    @staticmethod
    def shiftedTemperatures(T_units='°C'):
        """
        This function
        :param T_units: Units for temperature, default is '°C'
        :return:
        """
        header = ("Stream",
                  "Inlet temperature \n ("+T_units+")",
                  "Outlet temperature \n ("+T_units+")",
                  "Shifted inlet temperature \n ("+T_units+")",
                  "Shifted outlet temperature \n ("+T_units+")")
        colors = ("brown4", "midnight blue", "midnight blue",
                  "midnight blue", "midnight blue")
        for Header, Color, Cell in zip(header, colors, range(len(header))):
            tkinter.Label(shiftedTemperaturesTab,
                          text=Header,
                          background=Color,
                          foreground="white",
                          relief="raised",
                          borderwidth=3).grid(row=0, column=Cell + 1,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for i in range(len(T_in)):
            tkinter.Label(shiftedTemperaturesTab,
                          text=i + 1,
                          background="slate gray",
                          foreground="black",
                          relief="raised",
                          borderwidth=3).grid(row=i + 1, column=1,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for j in range(len(T_in)):
            tkinter.Label(shiftedTemperaturesTab,
                          text=T_in[j],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=j + 1, column=2,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for j in range(len(T_out)):
            tkinter.Label(shiftedTemperaturesTab,
                          text=T_out[j],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=j + 1, column=3,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for j in range(len(T_out)):
            tkinter.Label(shiftedTemperaturesTab,
                          text=Specifications[j][4],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=j + 1, column=4,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for j in range(len(T_out)):
            tkinter.Label(shiftedTemperaturesTab,
                          text=Specifications[j][5],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=j + 1, column=5,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')

    @staticmethod
    def cascadeTable(T_units='°C', Q_units='kW'):
        """
        This method
        :param T_units: Units for temperature, default is '°C'
        :param Q_units: Units for heat, default is 'kW'
        :return:
        """
        header = ("Temperature  ("+T_units+")",
                  "Heat  ("+Q_units+")")
        colors = ("midnight blue", "midnight blue")
        for Header, Color, Cell in zip(header, colors, range(len(header))):
            tkinter.Label(heatFlowCascadeTab,
                          text=Header,
                          background=Color,
                          foreground="white",
                          relief="raised",
                          borderwidth=3).grid(row=0, column=Cell + 1,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for j in range(len(temperatures)):
            tkinter.Label(heatFlowCascadeTab,
                          text=temperatures[j],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=j+1, column=1,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')
        for j in range(len(ptCascade)):
            tkinter.Label(heatFlowCascadeTab,
                          text=ptCascade[j],
                          background="khaki1",
                          foreground="black",
                          relief="ridge",
                          borderwidth=3).grid(row=j+1, column=2,
                                              columnspan=1,
                                              rowspan=1, sticky='nsew')

    @staticmethod
    def compositeCurve(Data='1', T_units='°C', Q_units='kW'):
        """
        This method
        :param Data: Name for the set of streams
        :param T_units: Units for temperature, default is '°C'
        :param Q_units: Units for heat, default is 'kW'
        :return:
        """
        global fig1, ax1
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        fig1.canvas.set_window_title('Hot and cold composite curves')
        ax1.plot(heatFlowSource, tempSource, 'r',
                 linewidth=2, label="Hot composite curve")
        ax1.plot(heatFlowSink, tempSink, 'b',
                 linewidth=2, label="Cold composite curve")
        ax1.axvline(x=min(heatFlowSource), color='k',
                    linestyle='--', linewidth=0.8)
        ax1.axvline(x=min(heatFlowSink), color='k',
                    linestyle='--', linewidth=0.8)
        ax1.axvline(x=max(heatFlowSource), color='k',
                    linestyle='--', linewidth=0.8)
        ax1.axvline(x=max(heatFlowSink), color='k',
                    linestyle='--', linewidth=0.8)
        ax1.minorticks_on()
        ax1.grid(which='BOTH', ls=':')
        ax1.grid(True)
        ax1.set_xlabel('Heat flow (' + Q_units + ')')
        ax1.set_ylabel('Actual temperature (' + T_units + ')')
        y_down, y_up = ax1.get_ylim()
        x1 = [max(heatFlowSource), max(heatFlowSink)]
        x2 = [min(heatFlowSource), min(heatFlowSink)]
        y1 = [y_up, y_up]
        ax1.plot(x1, y1, linestyle='-', color='orange', linewidth=2,
                 label="$Q_{Heating, min} = " + str(ptCascade[0]) + " \ " +
                       Q_units + "$", marker='o')
        ax1.plot(x2, y1, linestyle='-', color='green', linewidth=2,
                 label="$Q_{Cooling, min} = " + str(ptCascade[-1]) + " \ " +
                       Q_units + "$", marker='o')
        ax1.legend(loc="best", fancybox=True, shadow=True,
                   framealpha=0.98, frameon=True)
        fig1.tight_layout()
        return fig1.savefig('compositeCurve.jpg', bbox_inches="tight")

    @staticmethod
    def grandCompositeCurve(Data='1', T_units='°C', Q_units='kW'):
        """
        This method
        :param Data: Name for the set of streams
        :param T_units: Units for temperature, default is '°C'
        :param Q_units: Units for heat, default is 'kW'
        :return:
        """
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        fig2.canvas.set_window_title('Grand composite curve')
        ax2.plot(ptCascade, temperatures, 'r',
                 linewidth=1.2, label="_nolegend_")
        ax2.set_xlim([min(ptCascade), None])
        ax2.axvline(x=min(ptCascade), color='k', linestyle='-', linewidth=2)
        for k in range(len(temperatures)):
            x_1 = min(ptCascade)
            x_2 = ptCascade[k]
            y_1 = temperatures[k]
            y_2 = y_1
            x = [x_1, x_2]
            y = [y_1, y_2]
            ax2.plot(x, y, linestyle='--', color='k', linewidth=0.7,
                     label="_nolegend_")
            ax2.text(x_2, y_2, str(ptCascade[k]), verticalalignment='bottom',
                     horizontalalignment='left', fontsize=7, fontweight='bold',
                     label="_nolegend_", wrap=True, rotation=39)
            if ptCascade[k] == 0:
                ax2.axhline(y=temperatures[k], linestyle='--', color='blue',
                            linewidth=0.9, label="$T_{pinch} \ =$" + str(
                                temperatures[k]) + ' ' + T_units)
                ax2.legend(loc="best", fancybox=True, shadow=True)
            else:
                pass
        ax2.minorticks_on()
        ax2.grid(True, which='BOTH', ls=':')
        ax2.set_xlabel('Net heat flow (' + Q_units + ')')
        ax2.set_ylabel('Shifted temperature (' + T_units + ')')
        fig2.savefig('grandCompositeCurve.jpg', bbox_inches="tight")
        return fig2.savefig('grandCompositeCurve.jpg', bbox_inches="tight")

    @staticmethod
    def drawCascade(T_units='°C', Q_units='kW'):
        """
        This funtion
        :param T_units: Units for temperature, default is '°C'
        :param Q_units: Units for heat, default is 'kW'
        :return:
        """
        heatCascade = tkinter.Tk()
        heatCascade.title("Cascade of temperature intervsls,"
                          "energy balances and residuals")
        heatCascade.configure(background='white')
        emptyCellsTemperatures = len(temperatures)*[""]
        arraysTemperatures = list(zip(emptyCellsTemperatures, temperatures))
        cellsTemperatures = []
        for i in range(len(arraysTemperatures)):
            cellsTemperatures.extend(arraysTemperatures[i])

        emptyCellsCascade = len(ptCascade) * [""]
        arraysCascade = list(zip(ptCascade, emptyCellsCascade))
        cellsCascade = []
        for i in range(len(arraysCascade)):
            cellsCascade.extend(arraysCascade[i])

        for k in range(len(cellsTemperatures)):
            if cellsTemperatures[k] != "":
                tempLabel = tkinter.Text(heatCascade, foreground='black',
                                         background='white', relief='ridge',
                                         borderwidth=1, height=2, width=21)
                tempLabel.tag_configure("subscript", offset=-1)
                tempLabel.insert("insert", "T", "", str(k), "subscript",
                                 " = "+str(cellsTemperatures[k])+" "+T_units)
                tempLabel.configure(state="disabled")
                tempLabel.grid(row=k, column=0, columnspan=1,
                               rowspan=1, sticky='nsew')
            else:
                tempLabel = tkinter.Text(heatCascade, foreground='black',
                                         background='white', relief='ridge',
                                         borderwidth=1, height=2, width=21)
                tempLabel.tag_configure("subscript", offset=-1)
                tempLabel.insert("insert", "")
                tempLabel.configure(state="disabled")
                tempLabel.grid(row=k, column=0, columnspan=1,
                               rowspan=1, sticky='nsew')
        for k in range(len(cellsCascade)):
            if cellsTemperatures[k] == "":
                heatLabel = tkinter.Text(heatCascade, foreground='black',
                                         background='white', relief='ridge',
                                         borderwidth=1, width=21, height=2)
                heatLabel.tag_configure("subscript", offset=-1)
                heatLabel.insert("insert", "Q", "", str(k), "subscript",
                                 " = " + str(cellsCascade[k]) + " " + Q_units)
                heatLabel.configure(state="disabled")
                heatLabel.grid(row=k, column=1, columnspan=1,
                               rowspan=1, sticky='nsew')
            else:
                heatLabel = tkinter.Text(heatCascade, foreground='black',
                                         background='white', relief='ridge',
                                         width=21, height=2)
                heatLabel.insert("insert", "")
                heatLabel.grid(row=k, column=1, columnspan=1,
                               rowspan=1, sticky='nsew')
        heatCascade.mainloop()

    @staticmethod
    def initialGridDiagram(Data='1', T_units='°C'):
        """
        This method
        :param Data: Name for the set of streams
        :param T_units: Units for temperature, default is '°C'
        :return:
        """
        figGrid, axGrid = plt.subplots()
        figGrid.canvas.set_window_title('Initial grid diagram')
        axGrid.minorticks_on()
        axGrid.grid(True, which='BOTH', ls=':')
        for i in range(len(T_in)):
            axGrid.plot([T_in[i], T_out[i]], [i + 1, i + 1],
                        color='black', linestyle='-', linewidth=1.5)
            axGrid.arrow(T_in[i], i + 1, T_out[i] - T_in[i], 0,
                         head_width=0.1, head_length=1.3, color='black')
            axGrid.annotate(str(T_in[i]) + " " + T_units,
                            xy=(T_in[i], i + 0.92), ha='center')
            axGrid.annotate(str(T_out[i]) + " " + T_units,
                            xy=(T_out[i], i + 0.92), ha='center')
            axGrid.text(T_in[i], i + 1,
                        str(i), color='black',
                        bbox=dict(facecolor='white',
                                  edgecolor='black',
                                  boxstyle='square,pad=0.5'),
                        horizontalalignment='left',
                        verticalalignment='center')
        axGrid.set_xticklabels([])
        axGrid.set_yticklabels([])
        axGrid.invert_yaxis()
        figGrid.tight_layout()
        return figGrid.savefig('initialGridDiagram.jpg', bbox_inches="tight")

    @staticmethod
    def drawIntervals(Data='1', T_units='°C'):
        """
        This method
        :param Data: Name for the set of streams
        :param T_units: Units for temperature, default is '°C'
        :return:
        """
        figIntervals, axIntervals = plt.subplots()
        figIntervals.canvas.set_window_title('Intervals')
        axIntervals.minorticks_on()
        axIntervals.grid(True, which='BOTH', ls=':')
        for i in range(len(temperatures)-1):
            axIntervals.plot([temperatures[i], temperatures[i+1]],
                             [i + 1, i + 1],
                             linestyle='--', linewidth=1.5, marker='$▌$')
            axIntervals.annotate(str(temperatures[i]) + " " + T_units,
                                 xy=(temperatures[i], i + 0.8),
                                 ha='left', fontsize=7)
            axIntervals.annotate(str(temperatures[i+1]) + " " + T_units,
                                 xy=(temperatures[i+1], i + 0.8),
                                 ha='right', fontsize=7)
        for j in range(len(temperatures)):
            axIntervals.axvline(x=temperatures[j], color='black',
                                linestyle=':', linewidth=0.4)
        for k in range(len(T_in)):
            if T_in[k] > T_out[k]:
                colorLine = 'red'
            else:
                colorLine = 'blue'
            axIntervals.plot([Specifications[k][4], Specifications[k][5]],
                             [k + len(temperatures), k + len(temperatures)],
                             color=colorLine, linestyle='-', linewidth=1.5)
            axIntervals.annotate(str(Specifications[k][4]) + " " + T_units,
                                 xy=(Specifications[k][4],
                                     k + 0.98*len(temperatures)),
                                 ha='center', fontsize=7)
            axIntervals.annotate(str(Specifications[k][5]) + " " + T_units,
                                 xy=(Specifications[k][5],
                                     k + 0.98*len(temperatures)),
                                 ha='center', fontsize=7)
            axIntervals.text(T_in[k] + (T_out[k]-T_in[k])/2,
                             k + len(temperatures),
                             str(k+1), color='white',
                             bbox=dict(facecolor=colorLine,
                                       edgecolor=colorLine,
                                       boxstyle='square,pad=0.5'),
                             horizontalalignment='left',
                             verticalalignment='center',
                             fontsize=7)
        axIntervals.set_xticklabels([])
        axIntervals.set_yticklabels([])
        axIntervals.invert_yaxis()
        figIntervals.tight_layout()
        return figIntervals.savefig('drawIntervals.jpg', bbox_inches="tight")


# def file_save(event):
#     print('Click derecho')
#     f = tkinter.filedialog.asksaveasfile(mode='w', defaultextension=".png")
#     if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
#         return
#     Image.save(f)
#     f.close()

# Configuración de la ventana
WIN = tkinter.Tk()
#WIN.option_readfile('optionDB')
# WIN.iconbitmap(default='1.ico')
#WIN.wm_title("off")
WIN.title("HIntPy")
WIN.geometry('%dx%d+%d+%d' % (1268, 547, 0, 0))  # Ancho, Alto, X, Y
WIN.resizable(0, 0)
WIN.configure(background='lavender')

# Data Notebook creation
HInPyNotebook = ttk.Notebook(WIN, width=551, height=500)
# Tab for data
DataTab = tkinter.Frame(HInPyNotebook, background='gray25')
# Tab for stream specification
streamsSpecificationsTab = tkinter.Frame(HInPyNotebook, background='gray25')
# Tab for shifted temperatures
shiftedTemperaturesTab = tkinter.Frame(HInPyNotebook, background='gray25')

HInPyNotebook.add(DataTab, text='Data')
HInPyNotebook.add(streamsSpecificationsTab, text='Stream Specifications')
HInPyNotebook.add(shiftedTemperaturesTab, text='Shifted Temperatures')
HInPyNotebook.grid(row=0, column=0)

# Diagrams Notebook creation
DiagramsNotebook = ttk.Notebook(WIN, width=700, height=500)
# Tab for streams specification
streamsTab = tkinter.Frame(DiagramsNotebook, background='gray25')
# streamsTab.bind("<Button-3>", file_save)
# Tab for initial grid diagram
initialGridTab = tkinter.Frame(DiagramsNotebook, background='gray25')
# Tab for intervals of temperature
intervalsTab = tkinter.Frame(DiagramsNotebook, background='gray25')
# Tab for heat flow cascade
heatFlowCascadeTab = tkinter.Frame(DiagramsNotebook, background='gray25')
# Tab for cumulative enthalpies
cumulativeHTab = tkinter.Frame(DiagramsNotebook, background='gray25')
# Tab for gran composite curve
grandCompositeTab = tkinter.Frame(DiagramsNotebook, background='gray25')

DiagramsNotebook.add(streamsTab, text='Streams')
DiagramsNotebook.add(initialGridTab, text='Initial Grid Diagram')
DiagramsNotebook.add(intervalsTab, text='Intervals of temperature')
DiagramsNotebook.add(heatFlowCascadeTab, text='Heat Flow Cascade')
DiagramsNotebook.add(cumulativeHTab, text='Cumulative Enthalpies')
DiagramsNotebook.add(grandCompositeTab, text='Grand composite')
DiagramsNotebook.grid(row=0, column=1)


def LoadData():
    WIN.filename = tkinter.filedialog.askopenfilename(initialdir="/",
                                                      title="Select file",
                                                      filetypes=(("txt files", "*.txt"),
                                                                 ("all files", "*.*")))
    tkinter.Label(text=os.path.split(WIN.filename)[1]).grid(row=20, column=1)
    text = ScrolledText(DataTab)
    text.pack(expand=0, fill="both")
    progressBar = ttk.Progressbar(WIN, orient="horizontal", length=250,
                                  mode='determinate')
    progressBar.grid(row=20, column=0)
    progressBar.start()
    with open(WIN.filename, 'r') as f:
        text.insert(tkinter.INSERT, f.read())
    pinchStream(WIN.filename).streamData()
    pinchStream(WIN.filename).shiftedTemperatures()
    pinchStream(WIN.filename).drawStreams()
    pinchStream(WIN.filename).initialGridDiagram()
    pinchStream(WIN.filename).drawIntervals()
    pinchStream(WIN.filename).grandCompositeCurve()
    pinchStream(WIN.filename).compositeCurve()
    progressBar.stop()

def drawDiagrams():
    imageList = ('streamsCurve.jpg', 'initialGridDiagram.jpg',
                 'drawIntervals.jpg', 'compositeCurve.jpg',
                 'grandCompositeCurve.jpg')
    tabList = (streamsTab, initialGridTab, intervalsTab,
               cumulativeHTab, grandCompositeTab)
    for i in range(len(imageList)):
        # Imagen principal para modelo de absorbedor
        image = Image.open(imageList[i])
        image = image.resize((700, 500), Image.ANTIALIAS)
        photo = ImageTk.PhotoImage(image)
        label = tkinter.Label(tabList[i], image=photo, relief='ridge',
                              borderwidth=3)
        label.image = photo
        label.grid(row=0, column=0, rowspan=100)
    pinchStream(WIN.filename).cascadeTable()

def Info():
    """
    It show the information of application
    """
    Message(WIN, text='Hello',
            bg='royalblue', fg='ivory',
            relief='groove').grid(row=6, column=0,
                                  columnspan=7, rowspan=4, sticky='nsew')


def Quit():
    """
    It quits
    """
    imageList = ('streamsCurve.jpg', 'initialGridDiagram.jpg',
                 'drawIntervals.jpg', 'compositeCurve.jpg',
                 'grandCompositeCurve.jpg')
    for j in range(len(imageList)):
        os.remove(imageList[j])
    return WIN.destroy()


def About():
    """
    It shows
    """
    Pmw.aboutversion('1.0')
    Pmw.aboutcopyright('Luis Eduardo Correa Gallego \n 2019')
    Pmw.aboutcontact('Contact:\n' +
                     ' Universidad de Antioquia\n' +
                     ' Phone: (+57) 3012083600\n' +
                     ' email: luise.correa@udea.edu.co'
                     )
    about = Pmw.AboutDialog(WIN,
                            applicationname='HIntPy (Heat Integration on Python)')

def references():
    """
    It shows the references
    """
    sketch = """ References for coefficients in Antoine equation:

1. CRC Handbook of Chemistry and Physics, 75th 91th eds., CRC Press, Inc.: Boca Raton, FL: 19942011.

2. CRC Handbook of Thermophysical and Thermochemical Data, CRC Press: Boca Raton, FL: 1994.

3. Daubert, T. E. and R. P. Danner, Data Compilation Of Properties Of Pure Compounds, Parts 1, 2, 3, and 4, Supplements 1 and 2, DIPPR Project, AIChE: Newo Yrk, 19851994.

4. Encyclopedia Of Chemical Technology, 3rd and 4th eds., John Wiley & Sons: New York, 19782004.

5. Landolt, H. and R. Bornstein, Zahlenwerte und Funkionen ans Physik, Chemei, Astronomie und Technik, SpringerVerlag: Heidelberg, Germany, 19712005

6. Lange's Handbook of Chemistry, 13th, 14th, and 15th eds., McGrawHill, New York, 1985, 1992, 1999.

7. Lyman, W. J., W. F. Reehl, and D. H. Rosenblatt, Handbook of Chemical Property Estimation Methods, American Chemical Society: Washington, DC, 1990.

8. Nesmeyanov, A. N., Vapor Pressure of the Chemical Elements, Elsevier: New York, 1963.

9. Perry's Chemical Engineering Handbook, 5th, 6th, 7th, and 8th eds., McGrawHill, New York, 1973, 1984, 1997, 2008.

10. Poling, B. E., J. M. Prausnitz, and J. P. O'Connell, The Properties of Gases and Liquids, 3rd ed. (R. C. Reid and T. K. Sherwood), 4th ed., 5th ed., McGrawHill, Newo Yrk, 1977, 1987, 2001.

11. Samsonov, G. V., ed., Handbook Of The PhysicoChemical Properties Of The Elements, Plenum: Washington, DC, 1968.

12. Selected Values of Properties of Chemical Compounds, Thermodynamics Research Center TAMU: College Station, TX, 1997.

13. Selected Values of Properties of Hydrocarbons and Related Compounds, Thermodynamics Research Center: TAMU, College Station, TX, 1998.

14. Simmrock, K. H., R. Janowsky, and O. Ohnsorge, Critical Data of Pure Substances, Dechema Chemistry Data Series, Vol. II, Part 1 (AgC7) and Part 2 (C8Zr), Dechema: Frankfurt/Main, Germany, 1986.

15. Timmermans, J., PhysicoChemical Constants of Pure Organic Compounds, Vols. 1 and 2, Elsevier: New York, 1950, 1965.

16. Tsonopoulos, C., "Critical Constants of Normal Alkanes from Methane to Polyethylene", AIChEJournal, 33 (12), 20802083 (December, 1987).

17. Vargaftik, N. B., Tables on the Thermophysical Properties of Liquids and Gases, 2nd ed., English translation, Hemisphere Publishing Corporation: Newo Yrk, 1975, 1983.

18. Yaws, C. L., Chemical Properties Handbook, McGrawHill: New York, 1999.

19. Yaws, C. L., Matheson Gas Data Book, 7th ed., Matheson TriGas (Parsippany, NJ), McGrawHill: New York, 2001.

20. Yaws, C. L., Yaws Handbook of Physical Properties for Hydrocarbons and Chemicals, Gulf Publishing Company: Houston, 2005.

21. Yaws, C. L., Yaws Handbook of Vapor Pressure – Antoine Coefficients, Gulf Publishing Company: Houston, 2007.

22. Yaws, C. L., Thermophysical Properties of Chemicals and Hydrocarbons, William Andrew Inc.: Norwich, NY, 2008.

References for absorption calculation:

23. J. D. Seader, Ernest J. Henley, D. Keith Roper, Separation Process Principles: Chemical and Biochemical Operations, Thir Edition, John Wiley and Sons: 2011.

"""
    dialog = Pmw.TextDialog(WIN, scrolledtext_labelpos="n",
                            title="References",
                            defaultbutton=0,
                            label_text="References")
    dialog.insert('end', sketch)
    dialog.configure(text_state="disabled")
    dialog.activate()
    dialog.tkraise()


def manual():
    """
    It shows a brief tutorial for to use the application

    """
    sketch = """
DATOS DE ENTRADA
El usuario debe suministrar datos de porcentaje de humedad y de
densidad húmeda con base en cinco determinaciones de relaciones
densidad-humedad ejecutadas conforme a la normatividad. No se
permite el uso de un número menor o mayor de determinaciones.

UNIDADES Y CIFRAS SIGNIFICATIVAS
Sin tener en cuenta los métodos de cálculo utilizados para el
porcentaje de humedad y la densidad húmeda, se estipula el uso de
kg/m3 como unidades de densidad, tanto para el estado seco como
para el estado húmedo. La humedad debe expresarse con dos cifras
decimales, la densidad debe expresarse con dos cifras decimales.

DATOS DE SALIDA
Se genera un reporte mínimo en formatos variados, ilustrando cálculos
básicos para la densidad seca y calculando los valores de humedad
óptima y densidad máxima seca. Se genera también un gráfico para
la curva de compactación, donde se ilustran los puntos experimentales
y la curva de ajuste.
En la IGU (Interfaz Gráfica de Usuario) se generan elementos como:

- Casillas para los valores de humedad óptima y densidad máxima seca.
- Casillas para los parámetros de regresión, expresada como una curva
  polinómica de grado 2 (parabola).
- Gráfico de la curva de compactación, con su correspondiente etiqueta
  de creación.

Esta aplicación no utiliza atajos de teclado. Es una aplicación para
utilizar después de las determinaciones experimentales de humedad
y densidad húmeda, por lo que no admite datos de entrada como masas
de recipientes (vacíos, con muestra húmeda o muestra seca) o valores
para diámetros de moldes o sus masas respectivas.

"""
    dialog = Pmw.TextDialog(WIN, scrolledtext_labelpos="n",
                            title="MANUAL",
                            defaultbutton=0,
                            label_text="Manual")
    dialog.insert('end', sketch)
    dialog.configure(text_state="disabled")
    dialog.activate()
    dialog.tkraise()


# First step: create the menu bar
Barra = tkinter.Menu(WIN)
# Second step: create the menu's
Data = tkinter.Menu(Barra, tearoff=0)
Diagrams = tkinter.Menu(Barra, tearoff=0)
Help = tkinter.Menu(Barra, tearoff=0)
AvailableCompounds = tkinter.Menu(Barra, tearoff=0)

# Third step: create the commands in the menu's
# Commands for File menu
Data.add_separator()
Data.add_command(label="Load data",
                 command=LoadData, font=('Consolas', 10, 'bold'),
                 foreground='firebrick4',
                 background='lavender')
Data.add_command(label="Edit data", command=None,
                 font=('Consolas', 10, 'bold'), foreground='firebrick4',
                 background='lavender', state="disabled")
Data.add_separator()
Data.add_command(label="Quit", command=Quit, font=('verdana', 10, 'bold'),
                 background='red')

# Commands for Diagrams menu
Diagrams.add_separator()
Diagrams.add_command(label="Draw", command=drawDiagrams,
                     font=('Consolas', 10, 'bold'), foreground='firebrick4',
                     background='lavender')

# Commands for Help menu
Help.add_separator()
Help.add_command(label='Quick start', command=manual,
                 font=('Consolas', 10, 'bold'), foreground='firebrick4',
                 background='lavender', state='disabled')
Help.add_command(label='References', command=references,
                 font=('Consolas', 10, 'bold'), foreground='firebrick4',
                 background='lavender', state='disabled')
Help.add_command(label='About', command=About,
                 font=('Consolas', 10, 'bold'), foreground='firebrick4',
                 background='lavender')
Help.add_separator()

# Fourth step: Add menu's to the bar
Barra.add_cascade(label="File", menu=Data, font=('Consolas', 13, 'bold'))
Barra.add_cascade(label="Diagrams", menu=Diagrams, font=('Consolas', 13, 'bold'))
Barra.add_cascade(label='Help', menu=Help)

WIN.style = ttk.Style()
WIN.style.theme_use(
    'winnative')  # ('classic', 'default', 'clam', 'winnative', 'vista', 'xpnative', 'alt')

# Fifth step: menu bar will be in the bar
WIN.config(menu=Barra)

WIN.mainloop()

if __name__ == '__main__':
    pass