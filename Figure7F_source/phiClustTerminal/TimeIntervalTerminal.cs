using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;

namespace singleCellAnalisysTerminal
{
    class TimeIntervalTerminal
    {
        public delegate void UpdateProgress(object sender, EventArgs e);
        static UpdateProgress progress = null;
        //static System.Windows.Forms.Timer ti;
        static Timer ti;

        /*         public static void Start() { ti.Start(); }
                 public static void Stop(){ti.Stop();}

                 public static void InitTimer(UpdateProgress progress)
                 {
                    //ti = new System.Windows.Forms.Timer();
                    ti = new Timer();
                    ti.Elapsed += new ElapsedEventHandler(progress);
                    //ti.Tick += new EventHandler(progress);
                    ti.Interval = 3000; // in miliseconds
                 } */
        private static void RunEvent(object o)
        {
            if (progress != null)
                progress(null, null);
        }
        public static void InitTimer(UpdateProgress progressP)
        {
            progress = progressP;
        }
        public static void Start()
        {
            int timeout = Timeout.Infinite;
            int interval = 4000;

            TimerCallback callback = new TimerCallback(RunEvent);

            ti = new Timer(callback, null, timeout, interval);
            ti.Change(0, interval);
        }
        public static void Stop()
        {
            ti.Dispose();
        }

    }
}
