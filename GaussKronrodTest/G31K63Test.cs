using GaussKronrod;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using MultiPrecision;
using System.Collections.ObjectModel;

namespace GaussKronrodTest {
    public class G31K63Test<N> where N : struct, IConstant {
        static readonly ReadOnlyCollection<MultiPrecision<N>> x, w1, w2;

        static G31K63Test() {
            (MultiPrecision<Plus2<N>>[] x_plus1, MultiPrecision<Plus2<N>>[] w1_plus1, MultiPrecision<Plus2<N>>[] w2_plus1) = CoefGenaratorMP<Plus2<N>>.Coef(31);

            x = Array.AsReadOnly(x_plus1.Select(v => MultiPrecision<N>.Abs(v.Convert<N>())).ToArray());
            w1 = Array.AsReadOnly(w1_plus1.Select(v => v.Convert<N>()).ToArray());
            w2 = Array.AsReadOnly(w2_plus1.Select(v => v.Convert<N>()).ToArray());
        }

        public static void PolyTest() {
            static MultiPrecision<N> f(MultiPrecision<N> x) => 1 + x * (1 + x * (1 + x * (1 + x)));

            MultiPrecision<N> s1 = 0, s2 = 0;
            for (int i = 0; i < x.Count - 1; i++) {
                s1 += (f(-x[i]) + f(x[i])) * w1[i];
                s2 += (f(-x[i]) + f(x[i])) * w2[i];
            }

            s1 += f(-x[^1]) * w1[^1];
            s2 += f(-x[^1]) * w2[^1];

            Console.WriteLine("k_integrate");
            Console.WriteLine(s1);
            Console.WriteLine("g_integrate");
            Console.WriteLine(s2);

            MultiPrecision<N> expected = MultiPrecision<N>.Div(46, 15);

            Assert.IsTrue((s1 - expected).Exponent < -MultiPrecision<N>.Bits + 6);
            Assert.IsTrue((s2 - expected).Exponent < -MultiPrecision<N>.Bits + 6);
        }

        public static void SumTest() {
            MultiPrecision<N> s1 = w1.ToArray()[..^1].Sum() * 2 + w1[^1];
            MultiPrecision<N> s2 = w2.ToArray()[..^1].Sum() * 2 + w2[^1];

            Console.WriteLine("k_sum");
            Console.WriteLine(s1);
            Console.WriteLine("g_sum");
            Console.WriteLine(s2);

            Assert.IsTrue((s1 - 2).Exponent < -MultiPrecision<N>.Bits + 4);
            Assert.IsTrue((s2 - 2).Exponent < -MultiPrecision<N>.Bits + 4);
        }

        public static void WriteText() {
            using StreamWriter sw = new($"../../../../results/G31K63_n{MultiPrecision<N>.Length}.csv");

            sw.WriteLine("gauss/kronrod,x,w");
            for (int i = 1; i < x.Count; i += 2) {
                sw.WriteLine($"g,{x[i]},{w2[i]}");
            }
            for (int i = 0; i < x.Count; i++) {
                sw.WriteLine($"k,{x[i]},{w1[i]}");
            }
        }

        public static void WriteBinary() {
            using BinaryWriter sw = new(File.Open($"../../../../results/G31K63_n{MultiPrecision<N>.Length}.bin", FileMode.Create));

            for (int i = 0; i < x.Count; i++) {
                sw.Write(x[i]);
            }
            for (int i = 1; i < w2.Count; i += 2) {
                sw.Write(w2[i]);
            }
            for (int i = 0; i < w1.Count; i++) {
                sw.Write(w1[i]);
            }
        }
    }

    [TestClass]
    public class G31K63N8Test {
        [TestMethod]
        public void PolyTest() {
            G31K63Test<Pow2.N8>.PolyTest();
        }

        [TestMethod]
        public void SumTest() {
            G31K63Test<Pow2.N8>.SumTest();
        }

        [TestMethod]
        public void WriteText() {
            G31K63Test<Pow2.N8>.WriteText();
        }

        [TestMethod]
        public void WriteBinary() {
            G31K63Test<Pow2.N8>.WriteBinary();
        }
    }

    [TestClass]
    public class G31K63N16Test {
        [TestMethod]
        public void PolyTest() {
            G31K63Test<Pow2.N16>.PolyTest();
        }

        [TestMethod]
        public void SumTest() {
            G31K63Test<Pow2.N16>.SumTest();
        }
    }

    [TestClass]
    public class G31K63N32Test {
        [TestMethod]
        public void PolyTest() {
            G31K63Test<Pow2.N32>.PolyTest();
        }

        [TestMethod]
        public void SumTest() {
            G31K63Test<Pow2.N32>.SumTest();
        }
    }

    [TestClass]
    public class G31K63N64Test {
        [TestMethod]
        public void PolyTest() {
            G31K63Test<Pow2.N64>.PolyTest();
        }

        [TestMethod]
        public void SumTest() {
            G31K63Test<Pow2.N64>.SumTest();
        }

        [TestMethod]
        public void WriteBinary() {
            G31K63Test<Pow2.N64>.WriteBinary();
        }
    }
}