using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TwinCAT.Ads;
using System.IO;
using System.Diagnostics;

namespace ConsoleApplication1
{
    class Program
    {
        static void Main(string[] args)
        {
            TcAdsClient tcAds = new TcAdsClient();
            AdsStream ds = new AdsStream(16);
            bool scan = false;
            double input1=0;
            double input2 = 0;
            double freq = 0.2;
            double startFreq = 0.05;
            double endFreq = 0.2;
            double steps = 4;
            double periodNumber = 20;
            double time = 0;
            double amp = 2.5;
            double offset = 2.5;
            double deltaF = 1;
            double lastTime = 0;

            tcAds.Connect(851);

            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();
            if (scan == true) {
                freq = startFreq;
                deltaF = (endFreq - startFreq) / steps;
                while (freq<=endFreq)
                {
                    time = stopWatch.ElapsedMilliseconds / 1000f;
                    if ((time-lastTime) >= periodNumber/ freq)
                    {
                        lastTime = time;
                        freq += deltaF;
                        Console.WriteLine(time);
                        Console.WriteLine(freq);
                    }
                    else
                    {
                        ds = new AdsStream(8);
                        BinaryWriter bw = new BinaryWriter(ds);
                       
                        input1 = (amp * Math.Sin(2 * Math.PI * freq * time)+offset) * 3276.8;
                        input2 = (-amp/2 * Math.Sin(2 * Math.PI * freq * time) + offset) * 3276.8;
                        
                        bw.Write((int)input1);
                        bw.Write((int)input2);
                        tcAds.Write(0x4020, 0, ds);

                    }
                }
            }
            else
            {
                while (true)
                {
                    ds = new AdsStream(8);
                    BinaryWriter bw = new BinaryWriter(ds);

                    time = stopWatch.ElapsedMilliseconds / 1000f;
                    input1 = (amp * Math.Sin(2 * Math.PI * freq * time) + offset) * 3276.8;
                    input2 = (-amp / 2 * Math.Sin(2 * Math.PI * freq * time) + offset) * 3276.8;
                    input1 = 1;

                    bw.Write((int)input1);
                    bw.Write((int)input2);
                    tcAds.Write(0x4020, 0, ds);
                    Console.WriteLine(input1);
                }
            }


        }
    }
}
