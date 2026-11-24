"""Render all three figures from their companion CSVs (Python, NumPy, Matplotlib)."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

HERE = Path(__file__).resolve().parent
INK = "#172D35"
MUTED = "#475D67"
GRID = "#D9E1E3"
TEAL = "#086C70"
OCHRE = "#A34C28"
BLUE = "#365C86"
NEUTRAL = "#66737A"
WIDTH = 177.5 / 25.4

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.3, "text.color": INK,
    "axes.labelcolor": INK, "axes.edgecolor": "#93A3A8", "axes.linewidth": .65,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 7.7, "ytick.labelsize": 7.7,
    "xtick.major.width": .65, "ytick.major.width": .65, "xtick.major.size": 3, "ytick.major.size": 3,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": True
})


def data(n):
    with (HERE / f"Figure_{n}_source_data.csv").open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def column(rows, name):
    return np.array([float(row[name]) for row in rows])


def errors(rows, name):
    obs = column(rows, "observed")
    residual = column(rows, name) - obs
    assert len({r["gene"] for r in rows}) == len(obs) == 101
    return np.array([abs(residual.mean()), np.abs(residual).mean()])


def gain(rows, name):
    return errors(rows, "zero") - errors(rows, name)


def panel(fig, x, y, letter, heading, size=9.0, gap=.035):
    fig.text(x, y, letter, fontsize=11.2, weight="bold", va="top")
    fig.text(x + gap, y-.003, heading, fontsize=size, weight="bold", va="top")


def save(fig, n):
    fig.savefig(HERE / f"Figure_{n}.pdf", facecolor="white", metadata={
        "Title": f"Figure {n}", "Subject": "Forecasts of human periodontal ligament response",
        "CreationDate": None, "ModDate": None})
    fig.savefig(HERE / f"Figure_{n}.png", dpi=420, facecolor="white")
    plt.close(fig)


def fig1():
    rows = data(1)
    x, y = column(rows, "atlas"), column(rows, "observed")
    win = np.abs(x-y) < np.abs(y)
    assert win.sum() == 41 and (np.abs(x-y) > np.abs(y)).sum() == 60
    fig = plt.figure(figsize=(WIDTH, 3.10))
    panel(fig, .005, .977, "a", "Predicted and observed shifts")
    panel(fig, .513, .977, "b", "Error across the 101 genes")
    a = fig.add_axes([.092, .187, .361, .548])
    a.set(xlim=(-10.9, 11.1), ylim=(-26.8, 33.2), xticks=[-10,-5,0,5,10],
          yticks=[-20,-10,0,10,20,30], xlabel="Predicted rank change (pp)", ylabel="Observed rank change (pp)")
    a.axhline(0, lw=.75, color=GRID, zorder=0)
    a.axvline(0, lw=.75, color=GRID, zorder=0)
    a.plot([-11,11],[-11,11], color=MUTED, ls=(0,(4,3)), lw=.9, zorder=1)
    for mask, marker, color in [(~win, "^", OCHRE), (win, "o", TEAL)]:
        a.scatter(x[mask],y[mask],s=19,marker=marker,c=color,edgecolors="white",
                  linewidths=.35, alpha=.88, zorder=2)
    fig.legend(handles=[
        Line2D([],[],marker="o",ls="",color=TEAL,markersize=4.9,label="Closer than zero (41)"),
        Line2D([],[],marker="^",ls="",color=OCHRE,markersize=5.1,label="Farther than zero (60)"),
        Line2D([],[],ls=(0,(4,3)),lw=1,color=MUTED,label="Perfect agreement")
    ], loc="upper left", bbox_to_anchor=(.047,.917), frameon=False, fontsize=7.3, ncol=2,
       borderaxespad=0, handlelength=1.1, handletextpad=.4, columnspacing=.7, labelspacing=.35)
    methods = ["zero", "atlas", "atlas_shared_mean"]
    colours = [NEUTRAL, BLUE, "#657D94"]
    labels = ["Zero change", "Atlas", "Atlas mean only"]
    for i, (bottom, title, limit, ticks) in enumerate([
        (.534,"Absolute error of the panel mean",.86,[0,.3,.6]),
        (.082,"Mean absolute error across genes",8.55,[0,3,6])]):
        fig.text(.548, bottom+.324, title+" (pp)", fontsize=8.2)
        ax = fig.add_axes([.712, bottom, .268, .26])
        for yy, method, color, label in zip(range(2,-1,-1), methods,colours,labels):
            value = errors(rows, method)[i]
            ax.barh(yy,value,height=.42,color=color)
            ax.text(value+limit*.025,yy,f"{value:.3f}",ha="left",va="center",fontsize=7.9)
        ax.set(xlim=(0,limit),xticks=ticks,ylim=(-.55,2.56),yticks=[2,1,0],yticklabels=labels)
        ax.tick_params(axis="y",length=0,pad=4,labelsize=7.6)
        ax.spines["left"].set_visible(False)
        ax.grid(axis="x",lw=.55,color=GRID)
        ax.set_axisbelow(True)
    save(fig,1)


def fig2():
    rows = data(2)
    fig = plt.figure(figsize=(WIDTH, 3.14))
    panel(fig, .005, .974, "a", "Atlas transfer to stretch cultures")
    panel(fig, .548, .974, "b", "Prediction in a held-out donor")
    fig.text(.039,.848,"Day-3 tension atlas",color=MUTED,fontsize=8)
    fig.text(.582,.848,"Other four Hong Kong donors",color=MUTED,fontsize=8)
    a = fig.add_axes([.177,.267,.325,.464])
    xx0 = 0
    a.axvline(xx0,color=NEUTRAL,ls=(0,(3,3)),lw=.85)
    contrasts = ["Wuhan 24 h","Wuhan 48 h","GSE106887 12 h"]
    for yi, name in zip(range(2,-1,-1), contrasts):
        values = gain([r for r in rows if r["condition"] == name], "predicted")
        a.plot(values,[yi,yi],lw=1.3,color="#AAB9BF",zorder=1)
        a.scatter(values[0],yi,color=TEAL,marker="o",s=32,zorder=2)
        a.scatter(values[1],yi,color=BLUE,marker="s",s=27,zorder=2)
        for pos, xv in enumerate(values):
            a.annotate(f"{xv:+.3f}".replace("-","−"),(xv,yi),xytext=(0,6 if pos==0 else -7),
                       textcoords="offset points",fontsize=7.5,color=[TEAL,BLUE][pos],
                       ha="center",va="bottom" if pos==0 else "top",
                       bbox=dict(fc="white",ec="none",pad=.25))
    a.set(ylim=(-.57,2.58),xlim=(-2.24,.38),xticks=[-2,-1.5,-1,-.5,0],
          yticks=[2,1,0],yticklabels=["Wuhan 24 h","Wuhan 48 h","GSE106887 12 h"],xlabel="Gain over zero change (pp)")
    a.tick_params(axis="y",length=0,pad=5)
    a.spines["left"].set_visible(False)
    fig.legend(handles=[
        Line2D([],[],marker="o",ls="",color=TEAL,markersize=5,label="Panel mean"),
        Line2D([],[],marker="s",ls="",color=BLUE,markersize=5,label="Individual genes")
    ],loc="lower center",bbox_to_anchor=(.299,.035),frameon=False,ncol=2,fontsize=7.7,
       handlelength=1.0,handletextpad=.4,columnspacing=1.4)
    b = fig.add_axes([.646,.267,.313,.464])
    donor_scores = []
    b.axhline(0,ls=(0,(3,3)),lw=.8,color=NEUTRAL)
    for donor in range(1,6):
        scores = gain([r for r in rows if r["donor"] == str(donor)],"predicted")
        donor_scores.append(scores)
        b.plot([0,1],scores,lw=.9,color="#A4B9C0",marker="o",ms=3.6,markeredgewidth=0)
    means = np.mean(donor_scores,axis=0)
    b.plot([0,1],means,lw=1.8,color=INK,marker="D",ms=5.1)
    b.set(xlim=(-.16,1.16),ylim=(-.39,.64),xticks=[0,1],xticklabels=["Panel mean","Individual genes"],
          yticks=[-.3,0,.3,.6],ylabel="Gain over zero change (pp)")
    b.tick_params(axis="x",length=0,pad=5)
    fig.legend(handles=[
        Line2D([],[],marker="o",lw=.9,color="#93A9B1",markersize=3.5,label="Donor"),
        Line2D([],[],marker="D",lw=1.5,color=INK,markersize=4.5,label="Mean of five")
    ],loc="lower center",bbox_to_anchor=(.798,.035),frameon=False,ncol=2,fontsize=7.7,
       handlelength=1.35,handletextpad=.4,columnspacing=1.1)
    save(fig,2)


def fig3():
    rows = data(3)
    dest = ("Wuhan_24h","Wuhan_48h","GSE106887_12h")
    values = {c: [r for r in rows if r["condition"]==c] for c in dest}
    methods = ["zero","measured","measured_shared_mean","atlas","atlas_calibrated","atlas_shared_mean"]
    labels = ["Zero change","Hong Kong genes","Hong Kong mean only","Atlas genes","Calibrated atlas genes","Atlas mean only"]
    y = np.array([5.52,4.25,3.25,1.94,.94,-.06])
    fig = plt.figure(figsize=(WIDTH,3.18))
    panel(fig,.210,.978,"a","Error of the panel mean (pp)",gap=.031,size=8.8)
    panel(fig,.606,.978,"b","Mean error across genes (pp)",gap=.031,size=8.8)
    fill = {"ref":"#EEF1F2","less":"#D1E8E6","more":"#F2DFD1"}
    for i,left in enumerate([.245,.646]):
        ax = fig.add_axes([left,.175,.343,.657])
        ax.set(xlim=(-.50,2.50),ylim=(-.62,6.36))
        ax.set_axis_off()
        for j,title in enumerate(["Wuhan\n24 h","Wuhan\n48 h","GSE106887\n12 h"]):
            ax.text(j,6.02,title,ha="center",va="bottom",fontsize=7.7,linespacing=1.05,color=MUTED)
        for yy,method,label in zip(y,methods,labels):
            if i==0:
                ax.text(-.63,yy,label,va="center",ha="right",fontsize=8.0,clip_on=False)
            for j,c in enumerate(dest):
                er = errors(values[c],method)[i]
                diff = er-errors(values[c],"zero")[i]
                shade = fill["ref" if method=="zero" else "more" if diff>0 else "less"]
                ax.add_patch(Rectangle((j-.482,yy-.432),.964,.864,facecolor=shade,edgecolor="none"))
                ax.text(j,yy,f"{er:.3f}",va="center",ha="center",fontsize=8.4,color=INK)
    fig.legend(handles=[
        Patch(fc=fill["ref"],label="Zero-change reference"),
        Patch(fc=fill["less"],label="Lower error than zero"),
        Patch(fc=fill["more"],label="Higher error than zero")
    ],loc="lower center",bbox_to_anchor=(.5,.033),frameon=False,ncol=3,fontsize=7.5,
       handlelength=1.25,handletextpad=.45,columnspacing=1.3)
    save(fig,3)


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
