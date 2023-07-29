# list of produced files

# from subprocess import run, PIPE, Popen
import os
from time import ctime
import argparse
import re
import ROOT as r
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('var', type=str)
opt = parser.parse_args()

print(opt.var)

if opt.var == "Bpt":
    binsBP=7
    binsBs=4
    var_n = "pt"
    var_N = "PT"
elif opt.var == "By":
     binsBP=5
     binsBs=5
     var_n = "y"
     var_N = "Y"
elif opt.var == "nMult":
     binsBP=8
     binsBs=8
     var_n = "Mult"
     var_N = "Mult"
else:
    print("Invalid Input")

bp_pdf_list = [
    "henri2022/filesbp/signal_systematics_table_%s_ntKp.tex" % opt.var,
    "henri2022/filesbp/background_systematics_table_%s_ntKp.tex" % opt.var,
    # "BP/RawYieldFits/signal_systematics_table_Bpt_ntKp.tex",
    # "BP/RawYieldFits/background_systematics_table_Bpt_ntKp.tex",
    # "BP/RawYieldFits/general_systematics_table_Bpt_ntKp.tex",
    ]

bs_pdf_list = [
    "henri2022/filesbs/signal_systematics_table_%s_ntphi.tex" % opt.var,
    "henri2022/filesbs/background_systematics_table_%s_ntphi.tex" % opt.var,
    # "Bs/RawYieldFits/signal_systematics_table_Bpt_ntphi.tex",
    # "Bs/RawYieldFits/background_systematics_table_Bpt_ntphi.tex",
    # "BP/RawYieldFits/general_systematics_table_Bpt_ntKp.tex",
    ]

if not os.path.exists("syst_error"):
    os.makedirs("syst_error")

def get_pdf_syst(inFileList, outfile, hname, nbins, nshape):
    """
    get the biggest error from signal and background pdf variation
    return the root mean squared error
    """
    #2-d array with first line for signal pdf variation and second line for bckg pdf variation (biggest deviations only)
    error_sig_bkg = np.zeros([2, nbins])  
    
    for itype, file in enumerate(inFileList):
        with open(file) as fin:
            lines = fin.readlines()[2:]

        xgroups = re.findall(r'\d+\D+\d+', lines[0])
        xs = [np.mean([float(num) for num in re.findall(r'\d+', g)]) for g in xgroups]
        err = np.zeros([nshape[itype], nbins])

        for ishape in range(nshape[itype]):
            print( np.array(re.findall(r'\d+\.\d*', lines[ishape + 1]), 'd') )
            err[ishape] = np.array(re.findall(r'\d+\.\d*', lines[ishape + 1]), 'd')

        bigger_error = np.amax(err, axis=0)
        error_sig_bkg[itype] = bigger_error
    error_sum = np.sqrt(np.sum(error_sig_bkg ** 2, axis=0))
    # write output
    fout = r.TFile(outfile, "recreate")
    h = r.TGraph(len(xs), np.array(xs), error_sum)
    h.SetName(hname)
    h.GetYaxis().SetTitle("Total PDF variation error (%)")
    h.Write()
    print(np.array(h.GetY()))
    fout.Write()
    fout.Close()
    return

get_pdf_syst(bp_pdf_list, "syst_error/bp_pdf_%s.root" % var_n, "bp_error", binsBP, [3, 4])
get_pdf_syst(bs_pdf_list, "syst_error/bs_pdf_%s.root" % var_n, "bs_error", binsBs, [3, 3])

def make_line(item, array, suf = '', form = '{:.2f}'):
    entry = [item.ljust(22)] + [form.format(num) + suf for num in array]
    entry = [e.ljust(8) for e in entry]
    return ' & '.join(entry) + r'\\ ' + '\n'

def get_y(hist):
    return np.array([hist.GetBinContent(i)
                     for i in range(1, 1 + hist.GetNbinsX())])


def read_tracking_syst(inYield, nbins):
    """Return graphs and tables for eff, yield, and corrected yield"""
    fin = r.TFile(inYield)
    r.SetOwnership(fin, False)
    h_nom = fin.Get("hPtSigma")
    h_tight = fin.Get("hPtSigma_tight")
    h_tight_ratio = h_tight.Clone()
    h_tight_ratio.Divide(h_nom)
    g = r.TGraphAsymmErrors(h_tight_ratio)
    # g = r.TGraphAsymmErrors(h_nom, h_tight, 'pois')
    for ibin in range(g.GetN()):
        g.SetPointY(ibin, 100 * abs(g.GetPointY(ibin) - 1))
    # show percent error
    g.GetYaxis().SetTitle("Track selection error in %")
    ## syst table
    h_yield_nom = fin.Get("hPt")
    h_yield_tight = fin.Get("hPtTight")
    h_inveff_nom = fin.Get("hInvEff")
    h_inveff_tight = fin.Get("hInvEffTight")
    arr = np.array(g.GetX())
    error_tight = np.array(g.GetY())
    out_table = ''
    out_table += make_line('%s' % var_n, arr)
    out_table += make_line('raw yield', get_y(h_yield_nom))
    out_table += make_line('tight raw yield', get_y(h_yield_tight))
    out_table += make_line('loose raw yield', get_y(h_yield_nom))
    out_table += make_line('inverse eff', get_y(h_inveff_nom))
    out_table += make_line('tight inverse eff', get_y(h_inveff_tight))
    out_table += make_line('loose inverse eff', get_y(h_inveff_nom))
    out_table += make_line('corrected yield', get_y(h_nom), '', '{:.0f}')
    out_table += make_line('corrected tight yield', get_y(h_tight), '', '{:.0f}')
    out_table += make_line('corrected loose yield', get_y(h_nom), '', '{:.0f}')
    out_table += make_line('tight error', error_tight, '\%')
    out_table += make_line('loose error', np.zeros(len(arr)), '\%')
    print(out_table)
    return g, out_table

def get_tracking_syst(outfile, out_table):
    in_file_bp = "EffAna/BP/FinalFiles/BPPPCorrYield%s.root" % var_N
    in_file_bs = "EffAna/Bs/FinalFiles/BsPPCorrYield%s.root" % var_N
    g_bp, t_bp = read_tracking_syst(in_file_bp, binsBP)
    g_bp.SetName('bp_track_sel_error')
    g_bs, t_bs = read_tracking_syst(in_file_bs, binsBs)
    g_bs.SetName('bs_track_sel_error')
    fout = r.TFile(outfile, "recreate")
    g_bp.Write()
    g_bs.Write()
    fout.Close()
    with open(out_table, 'w') as fout:
        fout.write('\\PBp track selection systematics\n')
        fout.write(t_bp)
        fout.write('\n')
        fout.write('\\PBs track selection systematics\n')
        fout.write(t_bs)


get_tracking_syst("syst_error/syst_track_sel_%s.root" % var_n, "syst_error/syst_track_sel_%s.txt" % var_n)

def read_tracking_syst1D(inYield, nbins):
    """Return graphs and tables for eff, yield, and corrected yield"""
    fin = r.TFile(inYield)
    r.SetOwnership(fin, False)
    h_nom = fin.Get("CorrDiffHisBin")
    h_tight = fin.Get("CorrDiffHisBin_tight")
    h_tight_ratio = h_tight.Clone()
    h_tight_ratio.Divide(h_nom)
    g = r.TGraphAsymmErrors(h_tight_ratio)
    # g = r.TGraphAsymmErrors(h_nom, h_tight, 'pois')
    for ibin in range(g.GetN()):
        g.SetPointY(ibin, 100 * abs(g.GetPointY(ibin) - 1))
    # show percent error
    g.GetYaxis().SetTitle("Track selection error in %")
    ## syst table
    h_yield_nom = fin.Get("hPt")
    h_yield_tight = fin.Get("hPtTight")
    h_inveff_nom = fin.Get("InvEff1DHisvar")
    h_inveff_tight = fin.Get("InvEff1DHisTight")
    arr = np.array(g.GetX())
    error_tight = np.array(g.GetY())
    out_table = ''
    out_table += make_line('%s' % var_n, arr)
    out_table += make_line('raw yield', get_y(h_yield_nom))
    out_table += make_line('tight raw yield', get_y(h_yield_tight))
    out_table += make_line('loose raw yield', get_y(h_yield_nom))
    out_table += make_line('inverse eff', get_y(h_inveff_nom))
    out_table += make_line('tight inverse eff', get_y(h_inveff_tight))
    out_table += make_line('loose inverse eff', get_y(h_inveff_nom))
    out_table += make_line('corrected yield', get_y(h_nom), '', '{:.0f}')
    out_table += make_line('corrected tight yield', get_y(h_tight), '', '{:.0f}')
    out_table += make_line('corrected loose yield', get_y(h_nom), '', '{:.0f}')
    out_table += make_line('tight error', error_tight, '\%')
    out_table += make_line('loose error', np.zeros(len(arr)), '\%')
    print(out_table)
    return g, out_table

def get_tracking_syst1D(outfile, out_table):
    in_file_bp = "EffAna/BP/FinalFiles/BPPPCorrYield%s.root" % var_N
    in_file_bs = "EffAna/Bs/FinalFiles/BsPPCorrYield%s.root" % var_N
    g_bp, t_bp = read_tracking_syst1D(in_file_bp, binsBP)
    g_bp.SetName('bp_track_sel_error')
    g_bs, t_bs = read_tracking_syst1D(in_file_bs, binsBs)
    g_bs.SetName('bs_track_sel_error')
    fout = r.TFile(outfile, "recreate")
    g_bp.Write()
    g_bs.Write()
    fout.Close()
    with open(out_table, 'w') as fout:
        fout.write('\\PBp track selection systematics\n')
        fout.write(t_bp)
        fout.write('\n')
        fout.write('\\PBs track selection systematics\n')
        fout.write(t_bs)


get_tracking_syst1D("syst_error/syst_track_sel_%s_1D.root" % var_n, "syst_error/syst_track_sel_%s_1D.txt" % var_n)