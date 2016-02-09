import ROOT as root


def fit():
    gr=root.TGraph()
    for i in range(50):
        gr.SetPoint(i,i,i)
    
    func=root.TF1("func","[0]*x",0,100)
    can=root.TCanvas("fit","fit",800,600)
    can.cd()
    gr.Draw("PA")
    gr.Fit("func")
    

