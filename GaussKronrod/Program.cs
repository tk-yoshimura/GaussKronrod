using MultiPrecision;
using GaussKronrod; 

{
    (double[] x, double[] w1, double[] w2) = CoefGenarator.Coef(7, 1e-12);

    Console.WriteLine("x");
    foreach (double v in x) {
        Console.WriteLine(v);
    }

    Console.WriteLine("w1");
    foreach (double v in w1) {
        Console.WriteLine(v);
    }

    Console.WriteLine("w2");
    foreach (double v in w2) {
        Console.WriteLine(v);
    }
} 

{
    (MultiPrecision<Pow2.N8>[] x, MultiPrecision<Pow2.N8>[] w1, MultiPrecision<Pow2.N8>[] w2) = CoefGenaratorMP<Pow2.N8>.Coef(7);

    Console.WriteLine("x");
    foreach (MultiPrecision<Pow2.N8> v in x) {
        Console.WriteLine(v);
    }

    Console.WriteLine("w1");
    foreach (MultiPrecision<Pow2.N8> v in w1) {
        Console.WriteLine(v);
    }

    Console.WriteLine("w2");
    foreach (MultiPrecision<Pow2.N8> v in w2) {
        Console.WriteLine(v);
    }

    Console.WriteLine("w1_sum");
    Console.WriteLine(w1[..^1].Sum() * 2 + w1[^1]);

    Console.WriteLine("w2_sum");
    Console.WriteLine(w2[..^1].Sum() * 2 + w2[^1]);

    static MultiPrecision<Pow2.N8> f(MultiPrecision<Pow2.N8> x) => 1 + x * (1 + x * (1 + x * (1 + x)));

    MultiPrecision<Pow2.N8> s1 = 0, s2 = 0;
    for (int i = 0; i < x.Length - 1; i++) {
        s1 += (f(-x[i]) + f(x[i])) * w1[i];
        s2 += (f(-x[i]) + f(x[i])) * w2[i];
    }

    s1 += f(-x[^1]) * w1[^1];
    s2 += f(-x[^1]) * w2[^1];

    Console.WriteLine("integrate");
    Console.WriteLine(s1);
    Console.WriteLine(s2);
}

Console.WriteLine("END");
Console.Read();