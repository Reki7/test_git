using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace BioInfLib {
    public class BioInf_PA1 {
        private const string NucLetters = "ACGT";
        private const string ReserveNucLetters = "TGCA";

        public static List<int[]> GetCombinations(int k, int n) {
            if (k > n || n < 1)
                return null;
            List<int[]> combinations = new List<int[]>();
            int[] Narray = new int[n];
            for (int i = 0; i < n; i++)
                Narray[i] = i;
            int[] Karray = new int[k];
            for (int i = 0; i < k; i++)
                Karray[i] = Narray[i];
            bool finish = false;
            while (!finish) {
                //for (int m = 0; m < k; m++) {
                //    Console.Write(A[m]);
                //}
                //Console.WriteLine();
                var tmp = Karray.Clone();
                combinations.Add((int[])tmp);
                finish = true;
                for (int j = k - 1; j >= 0; j--) {
                    if (Karray[j] != Narray[n - k + j]) {
                        int val = Karray[j];
                        for (int l = j; l < k; l++)
                            Karray[l] = ++val;
                        finish = false;
                        break;
                    }
                }
            }
            return combinations;
        }


        public static List<int> FindPatternPositions(string pattern, string genome) {
            List<int> list = new List<int>();
            string search_string = genome;
            int cur_pos = search_string.IndexOf(pattern);
            int abs_pos = 0;
            while (cur_pos > -1) {
                list.Add(abs_pos + cur_pos);
                abs_pos += cur_pos + 1;
                search_string = search_string.Substring(cur_pos + 1);
                cur_pos = search_string.IndexOf(pattern);
            }
            return list.Count != 0 ? list : null;
        }

        public static string PrintPatternPositions(string pattern, string genome) {
            List<int> list = FindPatternPositions(pattern, genome);
            string result = "";
            if (list != null) {
                for (int i = 0; i < list.Count; i++)
                    result = result + list[i].ToString() + " ";
                return result.TrimEnd(' ');
            }
            return null;
        }


        public Dictionary<string, int> CountKMers(string genome, int k) {
            var dic = new Dictionary<string, int>();
            KMer kmer = new KMer(k);
            while (kmer.Next()) {
                List<int> list = FindPatternPositions(kmer.Text, genome);
                if (list != null) {
                    dic.Add(kmer.Text, list.Count);
                }
            }
            return dic.Count != 0 ? dic : null;
        }

        public string PrintFindMostFrequenlyKMers(string genome, int k) {
            string str;
            Dictionary<string, int> dic = CountKMers(genome, k);
            if (dic != null) {
                int max = dic.Values.Max();
                str = "";
                foreach (KeyValuePair<string, int> kv in dic)
                    if (kv.Value == max)
                        str = str + kv.Key + " ";
                return str.TrimEnd(' ');
            }
            return null;
        }


        public static string ReverseComplement(string pattern) {
            if(pattern == null)
                return null;
            string result = "";
            for (int i = pattern.Length - 1; i >= 0; i--)
                result += ReserveNucLetters[NucLetters.IndexOf(pattern[i])];
            return result;
        }


        /// <summary>
        /// Clump Finding Problem: Find patterns forming clumps in a string.
        /// </summary>
        /// <param name="genome"></param>
        /// <param name="k"></param>
        /// <param name="L"></param>
        /// <param name="t"></param>
        /// <returns>All distinct k-mers forming (L, t)-clumps in Genome</returns>
        public string PrintClunps(string genome, int k, int L, int t) {
            Dictionary<string, int> dic = CountKMers(genome, k);
            if (dic != null) {
                int max = dic.Values.Max();
                if (max < t)
                    return null;
                string str = "";
                foreach (KeyValuePair<string, int> kv in dic)
                    if (kv.Value == max) {
                        List<int> list = FindPatternPositions(kv.Key, genome);
                        for (int i = 0; i <= max-t; i++) {
                            if (list[i+t-1] - list[i] + k < L) {
                                str = str + kv.Key + " ";
                                break;
                            }
                        }
                    }
                return str.TrimEnd(' ');
            }
            return null;
        }


        /// <summary>
        /// Minimum Skew Problem: Find a position in a genome minimizing the skew.
        /// </summary>
        /// <param name="genome">A DNA string</param>
        /// <returns>All integer(s) i minimizing Skew(Prefixi (Text)) among all values of i (from 0 to |Genome|).</returns>
        public static string PrintMinimumSkewPositions(string genome) {
            List<int> list = new List<int>();
            int min = 0;
            int cur = 0;
            for (int i = 0; i < genome.Length; i++) {
                switch (genome[i]) {
                    case 'G':
                        cur++;
                        break;
                    case 'C':
                        cur--;
                        break;
                    default:
                        break;
                }
                if (cur == min) {
                    list.Add(i + 1);
                }
                if (cur < min) {
                    min = cur;
                    list.Clear();
                    list.Add(i + 1);
                }
            }
            string res = "";
            for (int i = 0; i < list.Count; i++)
                res = res + list[i].ToString() + " ";
            return res.TrimEnd(' ');
        }


        /// <summary>
        /// Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
        /// </summary>
        /// <param name="pattern"></param>
        /// <param name="text"></param>
        /// <param name="d"></param>
        /// <returns>All positions where Pattern appears in Text with at most d mismatches.</returns>
        public static List<int> FindApproximatePatternMatches(string pattern, string text, int d) {
            List<int> results = new List<int>();
            List<int[]> positions = GetCombinations(d, pattern.Length);
            foreach (int[] position in positions) {
                StringBuilder reg_pattern = new StringBuilder("(?=" + pattern + ")");
                foreach (int x in position) {
                    reg_pattern[x+3] = '.';
                }
                Match m = Regex.Match(text, reg_pattern.ToString());
                while (m.Success) {
                    if(!results.Contains(m.Index))
                        results.Add(m.Index);
                    m = m.NextMatch();
                }
            }
            results.Sort();
            return results;
        }


        /// <summary>
        /// Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.
        /// </summary>
        /// <param name="DNA">A DNA string</param>
        /// <param name="k">Length of k-mers</param>
        /// <param name="d">Mismatch level</param>
        /// <returns></returns>
        public List<string> FrequentWordsWithMismatchesAndReverseComplements(string genome, int k, int d) {
            var dic = new Dictionary<string, int>();
            //KMer kmer = new KMer(k);
            //while (kmer.Next()) {
            //    List<int> list = FindApproximatePatternMatches(kmer.Text, genome, d);
            //    int count = list != null ? list.Count : 0;
            //    List<int> list_rev = FindApproximatePatternMatches(ReverseComplement(kmer.Text), genome, d);
            //    count += (list_rev != null ? list_rev.Count : 0);
            //    if (count > 0) 
            //        dic.Add(kmer.Text, count);
            //}
            for (int i = 0; i < genome.Length - k; i++) {
                string text = genome.Substring(i, k);
                string text_rev = ReverseComplement(text);
                if (!dic.ContainsKey(text) && !dic.ContainsKey(text_rev)) {
                    List<int> list = FindApproximatePatternMatches(text, genome, d);
                    int count = list != null ? list.Count : 0;
                    List<int> list_rev = FindApproximatePatternMatches(text_rev, genome, d);
                    count += (list_rev != null ? list_rev.Count : 0);
                    if (count > 0)
                        dic.Add(text, count);
                }
            }
            List<string> results = new List<string>();
            if (dic.Count > 0) {
                int max = dic.Values.Max();
                foreach (KeyValuePair<string, int> kv in dic)
                    if (kv.Value == max) {
                        results.Add(kv.Key);
                        results.Add(ReverseComplement(kv.Key));
                    }
            }
            //results.Sort();
            return results;
        }
    }



    public class KMer {
        private const string NucLetters = "ACGT";
        private string _nucString;
        private int k;

        public KMer(int init_k) {
            k = init_k;
            Init();
        }

        public string Text {
            get {
                return _nucString;
            }
        }

        public void Init() {
            if (k > 15)
                return;
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < k; i++)
                builder.Append(NucLetters[0]);
            _nucString = builder.ToString();            
        }


        public bool Next() {
            return _next(_nucString.Length - 1);
        }
 
        private bool _next(int pos) {
            int newNuc = 0;
            int curNuc = NucLetters.IndexOf(_nucString[pos]);
            if (curNuc == 3) {
                if (pos == 0)
                    return false;
            }
            else {
                newNuc = curNuc + 1;
            }
            StringBuilder sb = new StringBuilder(_nucString);
            sb[pos] = NucLetters[newNuc];
            _nucString = sb.ToString();
            if (curNuc == 3)
                return _next(pos - 1);
            return true;
        }        
    }
 
}
