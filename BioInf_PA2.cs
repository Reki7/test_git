using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace BioInfLib {
    public class BioInf_PA2 {
        /// <summary>
        /// Protein Translation Problem: Translate an RNA string into an amino acid string.
        /// </summary>
        /// <param name="Pattern">An RNA string</param>
        /// <returns>The translation of Pattern into an amino acid string Peptide.</returns>
        public static string ProteinTranslation(string Pattern) {
            Protein p = new Protein();
            string peptide = p.RNATranslate(Pattern);
            return peptide;
        }

        struct RNAConstructorState {
            public int Position;
            public string AdvRNA;
            public bool IsAdvMatch;
            public bool IsRevMatch;

            public RNAConstructorState(int pos, string adv, bool is_adv, bool is_rev) {
                Position = pos;
                AdvRNA = adv;
                IsAdvMatch = is_adv;
                IsRevMatch = is_rev;
            }
        }

        /// <summary>
        /// Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
        /// </summary>
        /// <param name="DNA">A DNA string</param>
        /// <param name="Peptide">An amino acid string</param>
        /// <returns>All substrings of Text encoding Peptide (if any such substrings exist).</returns>
        public static string PeptideEncoding(string DNA, string Peptide) {
            List<string>[] encodes = new List<string>[Peptide.Length];
            Protein p = new Protein();
            for (int i = 0; i < Peptide.Length; i++) {
                encodes[i] = p.GetPeptideEncodes(Peptide[i].ToString());
            }

            string cur_adv, cur_rev;
            List<int> adv_match, rev_match;
            Dictionary<int, string> matches = new Dictionary<int, string>();

            RNAConstructorState cur_pos = new RNAConstructorState(0, "", true, true);
            Stack<RNAConstructorState> stack = new Stack<RNAConstructorState>();
            stack.Push(cur_pos);

            FileStream fs1 = File.OpenWrite(@".\PA2_2_debug.txt");
            TextWriter wr = new StreamWriter(fs1);
            foreach (List<string> e in encodes) {
                foreach(string s in e)
                    wr.Write(s + " ");
                wr.WriteLine();
            }
            do {
                if (cur_pos.Position < encodes[stack.Count - 1].Count) {     // position in the range -> check substring
                    cur_adv = stack.Peek().AdvRNA + Protein.UntranscribeRNA(encodes[stack.Count - 1][cur_pos.Position]);
                    //adv_match = cur_pos.IsAdvMatch ? BioInf_PA1.FindPatternPositions(cur_adv, DNA) : null;
                    adv_match = BioInf_PA1.FindPatternPositions(cur_adv, DNA);  //todo Fix the check (is previous IsAdvMatch spoiled?)
                    cur_rev = BioInf_PA1.ReverseComplement(cur_adv);
                    //rev_match = cur_pos.IsRevMatch ? BioInf_PA1.FindPatternPositions(cur_rev, DNA) : null;
                    rev_match = BioInf_PA1.FindPatternPositions(cur_rev, DNA);
                    wr.WriteLine("Adv[{1}.{0}]: {2} - {3}", stack.Count, cur_pos.Position, cur_adv, adv_match != null ? adv_match.Count : 0);
                    wr.WriteLine("Rev[{1}.{0}]: {2} - {3}", stack.Count, cur_pos.Position, cur_rev, rev_match != null ? rev_match.Count : 0);
                    if (adv_match != null || rev_match != null) {
                        // There are some matches
                        if (stack.Count < encodes.Length) {
                            // Not last level -> Go to next level (INTO)
                            cur_pos.AdvRNA = cur_adv;
                            cur_pos.IsAdvMatch = adv_match != null;
                            cur_pos.IsRevMatch = rev_match != null;
                            stack.Push(cur_pos);
                            cur_pos.Position = 0;
                        }
                        else {  // Last level -> Add the match and go to next substring (ADD + NEXT)
                            if (adv_match != null)
                                foreach(int pos in adv_match)
                                    matches.Add(pos, cur_adv); // Full string was matched
                            if (rev_match != null)
                                foreach(int pos in rev_match)
                                    matches.Add(pos, cur_rev);
                            cur_pos.Position++;
                        }
                    }
                    else {      // Go to next substring on the current level (NEXT)
                        cur_pos.Position++;
                    }
                }
                else {                  // x out of the range for current level
                    if (stack.Count > 1) {
                        // Not the first level -> Go to previous level (OUT)
                        cur_pos = stack.Pop();
                        cur_pos.Position++;
                    }
                    else {              // Last substring of last level -> EXIT
                        break;
                    }
                }
            } while (true);



            string str = "";
            foreach (KeyValuePair<int, string> kv in matches) {
                //str = str + kv.Value + " " + kv.Key + "\n";
                str = str + kv.Value + " ";
            }
            return str.TrimEnd(' ');
        }
    }


    class Protein {
        private Dictionary<string, string> RNA_codon;
        private List<AminoAcid> AminoAcidsTable;

        struct AminoAcid {
            public char _letter;
            private string Code;
            public string Name;
            public int Weight;

            public char Letter {
                get {
                    return _letter;
                }
            }
            public AminoAcid(string letter, int weight) {
                _letter = letter[0];
                Weight = weight;
                Name = _letter.ToString();
                Code = Name;
            }
        }

        public Protein() {
            Init();
        }

        private void Init() {
            RNA_codon = new Dictionary<string, string>();
            var codes = Regex.Split(Properties.Resources.RNA_codon_table, "\r\n|\r|\n");
            foreach (string code in codes) {
                var parts = code.Split(' ');
                if (parts.Length == 2) {
                    RNA_codon.Add(parts[0], parts[1]);
                }
                else {
                    RNA_codon.Add(parts[0], null);
                }
            }
            
            AminoAcidsTable = new List<AminoAcid>();
            AminoAcidsTable.Add(new AminoAcid("G", 57));
            AminoAcidsTable.Add(new AminoAcid("A", 71));
            AminoAcidsTable.Add(new AminoAcid("S", 87));
            AminoAcidsTable.Add(new AminoAcid("P", 97));
            AminoAcidsTable.Add(new AminoAcid("V", 99));
            AminoAcidsTable.Add(new AminoAcid("T", 101));
            AminoAcidsTable.Add(new AminoAcid("C", 103));
            AminoAcidsTable.Add(new AminoAcid("I", 113));
            AminoAcidsTable.Add(new AminoAcid("L", 113));
            AminoAcidsTable.Add(new AminoAcid("N", 114));
            AminoAcidsTable.Add(new AminoAcid("D", 115));
            AminoAcidsTable.Add(new AminoAcid("K", 128));
            AminoAcidsTable.Add(new AminoAcid("Q", 128));
            AminoAcidsTable.Add(new AminoAcid("E", 129));
            AminoAcidsTable.Add(new AminoAcid("M", 131));
            AminoAcidsTable.Add(new AminoAcid("H", 137));
            AminoAcidsTable.Add(new AminoAcid("F", 147));
            AminoAcidsTable.Add(new AminoAcid("R", 156));
            AminoAcidsTable.Add(new AminoAcid("Y", 163));
            AminoAcidsTable.Add(new AminoAcid("W", 186));
        }

        public string RNATranslate(string RNA) {
            string peptide = "";
            for (int i = 0; i < RNA.Length / 3; i++) {
                string str = RNA.Substring(i * 3, 3);
                if (RNA_codon[str] == null)
                    break;
                peptide += RNA_codon[str];
            }
            return peptide;
        }

        public List<string> GetPeptideEncodes(string peptide) {
            List<string> list = new List<string>();
            foreach (KeyValuePair<string, string> kv in RNA_codon) {
                if(kv.Value == peptide)
                    list.Add(kv.Key);
            }
            return list;
        }

        public static string TranscribeRNA(string DNA) {
            string RNA = DNA.Replace('T', 'U');
            return RNA;
        }

        public static string UntranscribeRNA(string RNA) {
            string DNA = RNA.Replace('U', 'T');
            return DNA;
        }


        public List<int> GetTheoreticalSpectrum(string Peptide) {
            List<int> spectrum = new List<int>();
            int petide_weight = 0;
            foreach (char c in Peptide) {
                petide_weight += AminoAcidsTable[0].Weight;
            }

            return spectrum;
        }
    }
}
