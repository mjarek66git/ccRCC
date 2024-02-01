using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;
using phiClustCore;

namespace Graph
{
    
    public delegate void ReferenceBoxChanged(bool change);
    public partial class DistanceControl : UserControl
    {
        public ReferenceBoxChanged refChanged = null;

        private bool HideSetup;
        public bool hideSetup
        {
            set
            {
                HideSetup = value;
            }
            get
            {
                return HideSetup;
            }
        }
        private bool HideRef;
        public bool hideReference
        {
            set
            {
                HideRef = value;
                if (HideRef)
                    HideReference();
            }
            get
            {
                return HideRef;
            }

        }
        public bool reference
        {
            get
            {
                return referenceBox.Checked;                
            }
            set
            {
                referenceBox.Checked = value;
            }
        }
        
        bool hideHamming;
        public bool HideHamming
        {
            set
            {
                hideHamming=value;
                radio1DJury.Enabled = !hideHamming;
            }
            get { return hideHamming; }
        }
        bool hideCosine;
        public bool HideCosine
        {
            set
            {
                hideCosine = value;
                radioCosine.Enabled = !hideCosine;
            }
            get
            {
                return hideCosine;
            }
        }
        public DistanceMeasures distDef 
        {
            get
            {
                if (radio1DJury.Checked)
                    return DistanceMeasures.HAMMING;
                else
                    if (radioCosine.Checked)
                    return DistanceMeasures.COSINE;
                else
                        if (radioPearson.Checked)
                    return DistanceMeasures.PEARSON;
                else
                            if (radioEucl.Checked)
                    return DistanceMeasures.EUCLIDIAN;
                else
                    if (fastHamming.Checked)
                    return DistanceMeasures.FASTHAMMING;

                return DistanceMeasures.HAMMING;
            }
            set
            {
                switch (value)
                {
                    case DistanceMeasures.HAMMING:
                    case DistanceMeasures.COSINE:
                        
                        if(value==DistanceMeasures.COSINE)
                            radioCosine.Checked = true;
                        else
                            radio1DJury.Checked = true;
                            if(value==DistanceMeasures.HAMMING)
                                referenceBox.Checked = true;
                        radioPearson.Checked = false;
                        label2.Visible = true;
                        break;
                    case DistanceMeasures.EUCLIDIAN:
                        radioEucl.Checked = true;
                        break;
                    case DistanceMeasures.FASTHAMMING:
                        fastHamming.Checked = true;
                        break;

                }

            }
        }
        public void FreezDist()
        {
            radioPearson.Enabled = false;
        }
        public void UnfreezDist()
        {
            radioPearson.Enabled = true;
        }        
        public DistanceControl()
        {
            InitializeComponent();
        }
        public void HideReference()
        {
            this.Size = new Size(this.Size.Width, 123) ;
            referenceBox.Visible = false;

        }

        private void radio1DJury_CheckedChanged(object sender, EventArgs e)
        {
            if (radio1DJury.Checked || radioCosine.Checked)
            {
                if (!hideSetup)
                {
                    label2.Visible = true;
                }

                if(radio1DJury.Checked)
                    referenceBox.Checked = true;
            }

        }

        private void radioRmsd_CheckedChanged(object sender, EventArgs e)
        {
            label2.Visible = false;

        }

        private void radioMaxSub_CheckedChanged(object sender, EventArgs e)
        {            
            label2.Visible = false;

        }

        private void referenceBox_CheckedChanged(object sender, EventArgs e)
        {
            if (refChanged != null)
                refChanged(referenceBox.Checked);

                       

        }



    
    }
}
