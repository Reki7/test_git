﻿using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using BioInfLib;

namespace BioInfLib_Test {
    [TestClass]
    public class BioInf_PA2_test {
        [TestMethod]
        public void ProteinTranslation_TestPattern1() {
            // arrange
            string pattern = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA";
            string expected = "MAMAPRTEINSTRING";

            // act
            string actual = BioInf_PA2.ProteinTranslation(pattern);

            // assert
            Assert.AreEqual(expected, actual);
        }

        [TestMethod]
        public void ProteinTranslation_TestPattern2() {
            // arrange
            string pattern = "AUGCCCAUGGGAUUAGUGUGGCACAAACAAGGACCACUAGAAAGGAUAUCUAUAAGAGGAGUAAUAGGAGUUAGGAGCGGGUAUAACGAAACCAUUCGAAGGAAUUGGGUCAUGUUAGUAAGUAAAAGCGCCUUAUUCGUAUCCACAUGCUGCCAAUGUAACCCCCCUUACCUGACUUGUUAUAAGCAGUUGAAGAGUCCAGACGUGACACGUUUUGCGCGCGCUCAUGACAUGGAUCAUUUUAGAGACCACACUCAUAUGGCUGGGAGGACCAACUUGGAACAAACGUUUUGUGCUCAACCCGUACAUCUGACUAUGGACCUAGAGUAUUAUCAGGACCCUCCAGUGGCUUAUGUACUACAGUACAUGGUGCCGAGGCGAAUGCCAGCUCUGUCAGUAAUCACGAAUCCUCCCAACCAAGAAUUGCAUUCACUCUGGUCAUAUCUGCAUUUGUCGGUUAGUCAAAUCUCUCGGGGGCGUCUAACUUGCGUGACGCUGCUACUAGUGAUCUGGGACAACCAUUUAUUACGGUACUCCUAUCUUUGCCAGGCCUUCGAAAUGGUUAUGAGCUCGAACAUCCGCUCCUUGGGAGCUGAGACGUUAAAGGUGUUAACAGAUCCGAAUAGAGAGUGCGCUUCGGUCACCAUCCCACAUCUCACUCAGCAAUCCGCUCUGAGUAUUAGAGCACCCCCCAGGUACCUUCACACUGAAUUACGGCGUUGUAGAACAUCUUUAGCACGGACGACCCUUCGGUUCGUGUCCAGUUGUACCUCUCGAUUUGGUGCUUCGAUUCCGAGCUUUCAUGAACUCUACUGCUCGCUAUACACAGGGCUCUCCUUCGUUCAUAUAUCUCGAUCAAUUAAGAACGCGGUUGAUCAAACAUGCGGCAAACCUAUCGACAAUACUCAGAUUGUGAUAAGCUUCCAGAGCUCUCGUCUCCAUCGAUGCAGUCUCAUAGCUCGUCUACAGGCCGACCGUACUUUCUUCCUUUAUAAGCGAGGGUGGUGCCGGGGCGUGUGCUCUAUAGUUAUAUACUUGACCGGGAUCUGCUACGUCAAACUGUGUUGUUACGCUACGGGCAUUAGCACCGUUAGCGACAGCUAUAAACCAAGGGGUUACGUUCAUCUUCAGGUAUUUCUCUUCCGACACGUGACUGUGUACGUAAUUCGUACAUAUGAGAAUUCAGUUACCGGGGCUUCACGAUAUCAAGCGGACUACAUCAUGGGGCGCGACUUACCUACUCUCGUAGAUAGUUUGCCAGUGGUGUCGUAUAGUUGGAACCAGUACCUGCUUCUAGUCUGUGGACCACGUCGGGGAGACCUCUUAUGGCCCCUAUUAGGCAGGGAUCCAACAGUUAAAUUUAAACCGCGGACCUACCGUCAACGGAGGAGUUUUACAGCUUUCAUGACGCGGAAGUGCGUUUCACUAUGGGCGCAAUUGACAGCUCGCCGUCGGACAAGCCGCAAAUUAAUUCGGGCACUAACUCACGGGUGGAAAAUGCGAAUGUCAAGACACGCAGUGGUCUAUACAGGCAUAGAGGGGAUCCUCGGGUCGAAAACAGAGUUAACGGUCCCAACUCAUCCGCUGUACAUCUUGCUCCCGUGCUCGGGUCCAAGUUGCGUGCUCAUGACUGGGGUCCGACCAUUUUGGUCUUCGUCCCAAGAUAGUCGUUGGUCAACAACCAUCAGGCCGGAGGGUGGUGUACUGAGCGUAAGUUGCAGACAGAUAGCCAAUAGCACCUUCUCGGUGCUCGAGUCAUUAGGUCUCUUUAUAGAAGUUCGGCACGGCCACGGGAAAAUCCCGCUCUACCGUAGCUCCACCUGUAGUAAUUGUUCACAUGUCUGUCAGUCUAAUGAGUGGACAGCGUGGUUCUUGAACUCCCCUGCGGCCGGGCCGAACCAGUGCCAAAUUGUUUAUAAUACGAAAUACUGUAUUGCAGGGUACGCCCCGAGCCCGUUACUAAGCACAGCCCGCGCAGCGAGUUAUCGAUAUAAAUGCUCGUGGUACUUCCUUCUUCUCUUUUUGACGUGUUCCCUCAUCUCGGGGCAAUUCUUGGAACCGGGGAUUAUGAAGUGUGCACUGGGCAUGCUAAAGGUAGGCUCGGAAUGCAAUCUCGAUGAACCGCACGUUAGUUCCACCGCGGGUUACCACGGCCACGUGAACACCCUCUGUUGCAUUAUUGCAGUUAGCCCGCCGCUAAGCCUGCGGGUUAAAAGACCGGCGGGGAGAGUUGUCUUCGAAAUAGGUCCCGCAUGUGAUCAGGACGGUGGUCAACUUGCGUCAAAUAUAGUGAUAACUCGUAAUCCGACUGGGUGGCACCCAUGGACUUCCCGAGUCAGCGCCCUUACCAUACGCUCACUACAAUGCCCACAGUGUGUGAACGAAGGCUGGCACUGGUCGCACCCUGAGCCAGGACAUAUGGGACGAUACUUGCACCCAAACUCAAUCCCUUACGAAGCCUGCGUCUCUACCUAUACAAUUGCCGAAAAAUUCUCGGUUCAUUUUAGUACCCGGAGCAUUGACUGGGAGUGCCCGGUCCGAAUACUCGCGACGUACACAAGACAGCUGGAGUUGCGUUGCCUUACCGGAACGUGCUUUUCCUACUGGAAUACCCCAAAGAUUCCUACGGUAACACCCCAUGUAGCACGGAAAUCUAAGGAGACUGACCGAAACGCCUCGAUUCAGCCCCCGAUGGCUGCUACGGUACCAGCAACAGGCUUAUCGAAGCAAUGGAUUGUCAUCCAAUCUUCCAGGGAUCUCCGGACUGACGAUGAGAUGGAUACCACCUUUCAGACACAUCGGUAUGGCUUCGGAAUAUCAACAACGGAGUGCUACACUCUACCGAGCCUAUACUCUGCCCCAAAGAGUAGAAAUACCCACACUUGUCACGCAUCGCCUUCCAACGAUAGUUGGUCCUGCUUCAAGCUCUGUGGAGUAAGUAUACGUGAUAGGAAAAGACGAGUCUGGACUGCUUCGCGAUGUAAUGACCUAGCGUCAGAAAAGGGUGGCCGGCAUGGUACACCACCCACAAGAAUCCCACCGUCGUUGCCUAGAAAAACGACACAGACUGGUGCAGCGAACCUAGCUUACUGGUCCCAGUUCAUUAACGGCGACAAAUCGCCGGAUGGCUAUACUAAGCAUGUUGGCACGGACGCCAUAUAUGCAGCGCGUGUUCAAGCGCACAGCCUGACGUGCUCCCGGCUAUCCAGAAUUGUGUGUAGUGGAUCCGCCGUGGGUUGGGAUGGAGACGAUUUUCGUGAGUCGGACGUUGACCGAGAAUGCGAGGUUACGAGUGGCCCAAGGCCCCAGAGAUUUUGUAUUUACAUCAUGUCCUCCGUUUCGACCGACGGUGGAGUAUGUCGAAAGAUCUGGAUAUUCCACCUCUUAUCCAAUUCAUGUAGGAGCUUGUUUUCCCCGUCAUCAGAUAGGGAACAGAAGCUCGGAAUAAUGGCGCAGGCAGAGGUGUUCUGUCACGCCCUCUCAAAAUUUCUGUUUCGGUUGUGUUCUAUUUCAUCGCUCAGGGCGCUCAGUUCUUAUCCCAGAUCGGGGCGAGUUCUCUCUAUUGUGACAAACAUACCGCAGUGGGACAGCCGCCCUAGGCUACGAAUUGACUACCAUCAAUUUACACCAAUUUUUCAUAACCCGCCAGGUAAACUCCUGCAGCCUAACGAUACGUGCCUGGGCUGUUCAAACUGCCUCCAGAGCUCGCCUCUGCUUUUUACCCGGGCAGGUCUUUUGGGCCUAUUCAGAACCAUGGCCCCUAUAACAUACGGCUGCAUGACUCGCUCUGCUUGGUCUUCGUUACAGCGUGGAGGGCCACUCGACGUUAUCCUCUCUACGUCGAGAGCUACCCUAUUCAAGCGGUGUGUAAGAUCACGAUUCAGAAAUCCAUGUGGUAUUAAUGUUCUCACAAACUUUCACGAACCGCGAGCACUUCUAUUUAGCGCUCCACGCGGGAAUCCGAGCUAUACGGUCAUGUUCGCCAGUACGGCAGCCACGUGCGAGCGCCUGGCGUUAGCACGUCAUCGUACUUAUUCCCGAUUCCUCUUUGAACGUACUAGAAAACAUGCGCUGAAAGGGAGAACAAGUGAUCAAGACGGCCUAAGUGCCGGUGUCACGUGGAAGGGCGUGAAUUCUGGCCGGUACCUUUAUAGCACGUGUCGUCAUGGGCUGCGCAGGUUGAUCCGGACGAUGAUGAAGCAAAUUGGAAGACUAUUUUAUAGAACGAGGUGCACGCUUCUCGUCCACAGAGACCCUGGCGCGGAACACUACCUUAGUACGCCAAUAGUAAAGCAAGGACGCGGCUACCACCGGCGUACACGACUUCCAUUUCUGCUUUUUUACGACCAAUGCCCAGCGCGCGGAGUGCCGCUAUUAGUACGCUUAGUGACCAAUGCUUUAUCGCCUUAUAUUACCAUGUCAUAUGCAUGUCAGACUUGCAACAGAGCCCAGGAGAGCGUCAGCCUAUUACCGAGUACAUGUGCUUACCCCUACCGGUCAAUCAAUUCCGAUCCUGUUGGAGCUUGGAAAGUGGGCUCCCUAAAUACGCACUACAUUUGUAAGUUACCCAAGGUAUCACCGCCGAGAUGGGCAAGCUCACAGCAUAUCGCGAGCACAGUCGUUAGAAGCGAUACCCAUAGUAGAGACGAGUCGAUCAAUGGACAUGGCAGGCUUAAAGUUGCAAUCAUUGUUUGUAGGGUAACCGCAAGAAGUGUAGCAGAAGUCGGAGGUUACUGGAAAAAAAGACACAGGAGUCCACGAAGUAUAUUUGAUAUGAAAAAUAGAGACACUUGGCACGAUGGUGGGCGCGGCCCUACAUCAUCGUCGGACCAUUAUUAUCCCAUAACAGGCAUUCCUGUUCUCAGUACUCCCUCCCGGAGGCCUCAAGCUCGGCUCACCAUUAAAUCGUCAUCAAUUUGGGCUUAUCGAGAUCAACCGAUGCUCACACAGCCAUGCAUUACGUACCCGGUACUGCUCUUGCAAUGGGCCGUAGCUUCGCGAAUGAGUGCAAUAUCUGGCAUUCGGCUUGCGGGUCUGUGGGGAGCAGCUAAUGCACCAGUCCAACUCUUCCUCUACCUCAGUCAGUGCGCACACUUAGGAACUUAUUGUGUAUAUCACGAAAUGGGCCCAAUGUAUAUGGUUCUACCAAACCAACGGCGACUGCGACGCUUUAGUCUUGGAGGACAGAUUCUAUGUUACCCCUUCAGGACGACUUACAGUCGGGCUGUCGUUAGGCCUUAUCCUGAUACGCCGCAUGGUUACAGCGGUAACGGCGUAAACACCGGCUUCCUAAGACGCUGGCCUGGACCUGUCCCGUGUCUGUGGGUUAAGGGUGUUGCCCACCACGCGUCUAAGCUACGGUGGAAGACACUUGCUUGUACUGCGUUAUCCCUGCCUGCCCUACCUAGUCACCCUGGAUUAUUGAAAAGCCUUCUCAAAAGAGGCGAAGUGAAAAAACUGGAGGCAGCUCGCACCGUCUUGCGAGCUUUACAACCGCACCGAACAAAUGGCGGAUAUACUAAGGAACUCGGUCAGGAGGAGGAGCACUGGUUCAGGUACUCAUUCGGAGGUACGAUCGGGGCCUUGUACGAUCGUAAGGUUGCUGUGUGGUACCGCAGCUCACAGGUAAUCUACAAUGGUCCUAUUUAUCGUCUAGUACAGGAGCCGAGUUGCGCCUCUCCGAUCACCUUAUUCAGCCGAGUACUACUUAAUGUUGUAUCUUAUCUAGAGAUAACGGUUAAACGACAUCUGUUCACACCAAACACCAUCCCCGGACGAGGCCUUUCUCAUCGAAAGGGGGCAUUAUUUUUCCUAAUACCACCCUAUCCAAGCUCGUUCGAAUUAAUAAGAAAACCCUACCGCCGCACACAGCCUGCUUUUCCGGGAUCCCUAGGCGAUGAUCCCUUGCUCCCCUGUAUCUGUCGGUCUAAACUGCGACUCAGACGGAUAACUGUGUAUGGUUCAAGCCGGCCAGUUGAAGACGUGGGAACCACGGAGUCUAGAUGUCCAAAAACAUUUUACUCUCUCCCGCUCGCGCUGUGUUGGUGCCGGAGUACGGGGGUGGGCAUCUACGGCUAUCUGGAACCGGACUACCCGUAUUAUACGAGGCGCGCCCCCGCCAUACCAAUUUGGAGACACCCUGUAUAUCGCGCUUUCGUACAGUUAACGUUAGCGCAGGCGACGGUACCCGGGCCGAAGCCGUUUUUAUUGCUACCCGUGGGCACCGGACCGUUUAUAAAAAUUUUCGGCCGCCACGGUUUUGCAGUUUUAGCUGCCGUACGAUAUUGUACGCUGUGCUCGGACCCUAUACGGGGACCCGACAGCCGACUUAGCCUAGUAACAGAGAGUGCUCGCCAGGUUCUCACUCCGGUCGAGUAUACCGAUCUGCAACUGUGUUCACUGGGGAGCAUUCAAAAGUCAAAACGGCUCUCAAGACUAUCCCUAACCAAAUUUGAUGAAUUAGCCGCCCGAACAAAUGCUCGGUUCGUGCGUCCCGUAGGAUAUCAACAGGUUAAGCCCCGGGUGAUUGAAACUGAAACUCCUUAUAGGACGUAUCGGCCAAACCCGCCUACAGAUGUAGACGAGGAGCACAAACAUCGUUUACUAGCUGUCAAAAGAUACAUAGGUUUCAGAGGGACACUGUCGGACAGACCCCAGGCGGCGCUCAUUAAAGGUGAAUGCAAGGCCAGCAGGUCGGGGCAUCUUCUUGAAAUUCAUAGAGGAUCAGAGAGACUUGAUUGCCGUAUCACGCUAACGCUACUUAAGGACCCACUUGGCUACUUUUACGGAGUGGACAGUCAGCCACCGGCGGUGCAUUUGAGUUACGACGCGUUCUCGCAACUCAUCAAAAACAGAUGUUUUAGAAAUUACAAGUCGCGAGCUUUCCAAGGCACAAAUGUAGCACGAAAAUCAAUAUGGGAUCGGUCAAGUGACCCUGGUUUGAGGGAGCCAUUGUGUUGGACACUAGACAAGCGCUGGGCCAAGUACACACAUCAACCUCAACCCUGCCAGCCGGUUCCCACCUCCAAUGUAUUUGAAUGCAGGCGUGGCCAAACCGAGGUACGACAUAACGCGUUGGUCUAUAGUCCAACCUUCGCGAAUUUUCGUGUAUCAGCGGACAGUACUGUCUUACUUCCAGCCACCGGGUCAUUAAUCGUUCCAUCUGGAAGGUCCGAGCGUAACCCGGAUUAUCACGCGCCGUGCUCUUUAUGUCCCAACCUCAACUCUGGUUUGCCGGGUAGGAUAGUAGGGACGCCUCGAACCGACAAGCGGCGGGCCGUCACGCAACUCAGGGUGUCCGAAGCAAGUUGCCCUCAGCGGAAUAGAUCGCCUGACUCGUGCGCACUAGGUAAAUGCCCGCUUCCGUUAGAUAACUCUACAAGCCGAUUUUGCGGCCCUUGCGCACAUUAUUCAGACCUAACCACAUUUUUACGGCUCUUAUUCUCACAGACUUCGCCCGAUGGCAGGUUUUUCCAUCAGAUAAUAGCGGUACUCAAAGUGGGUAAAUGCAUUAGAGCUACUGCAGCAUUGAACCCCGCCUCGUACCUUGCCAGCCGUGCUGCACAGUCUGGAAAUAUUGUGGAAGGUCAUCGGGGACCGAUUUUAAGCAGGAGCAUGCUAGUAAAGCGAAACUGGCGUUGUCACGGGUCCUUACCAUGCAAUAGGACGACGGGAGAAACACUUACGAUAUUAGGCAUAAUGUACGGAUCCGUGGGUGGCUCCAGAAUAUUACAAGCAGUAGUGCAUUCGGUGGCUAGACGCGCCCCGCAGGCGACCGAUGGUUUUGGGCAUUUGCUGACCAAACCGUUCGUCAGAGCCAGAACAUCAGAAGAGGAUAAAGGGUGGGUGGUGACUUCUAUAUUCCCUCGACGGAACGCCUAUGACCCACUCGGUAAAGUGGGCCCAACAUGGCGUACGCCCGGCUAUCGCAUGAUACCUCAAGCAGAUAAGUGGGGGAUCGGACCGACUAACGGUGAUACAGCCAACUGCGUAGGCGGGUGCCUCAGGGUAAUUUGCCCGUCCAGGGCAUCGAAGCAUCACGCAGCAAUAUGCUUAUCCGGCACAAGCUUCGCACGUGUGGGAAAAGGAAGAGAGGAAAAAGUCCGUCAUACGGAGUCCCAAGCGUGGACACGCUUUAAUCUGGGUAAUGCGAGGCGAGGUGGGGUGGACCAGAUUUUACGGAAAUUAACCAUCGAUAGGGGUGGGCUACGUGGGCACGCGAGGAUAAGUACGAAAUACCGGCCCCGCCUACAGAGGAAUGCCACAAUGGCCGGGAAGGAUACCCGACAUUAUUUGAGUGCUUUCAAGUUCGAACUUGCCACGCGAGCGUACAGAUCGGCGCUAAUUCUCUGGGAACUUAUCCAGAGGAGUACUAGAACGCUGGCUUAUUCUACGGGGCGGAUUGAACUCCCAACAACAGUGCCCCACUCUUCGCCCGGUGUUGUGCUAUUAGGAGUAGUGAAUUUAAGAGGAGAAGCCACACAGCAUGCCCUAUUCUCCUAUCUGGAUUCUACUACGCCCUAUGAGUGCUUAGACAUAGCCCCCAAUAGCUGGUUCACGUCCGAGUCAUCCUGGCCACGUAGCGAUGCAAUACUAUCAUACCGGUUGGUUAAUGACCCGGCUGAACUUCUUAGUGCAGGGCCAUCAACGCUAUCACACCCGCUCUGCGGAUCCAUCGGCCGCUCCUUAACACGAGGAGACAGAUUAACGAAGGACUACUACGUAGUACACAAUUUGAGCUCCUCGCUUACGCGAUUGGACCGUCAUAGUUCAGCGACGCGCAGCUUGCGACCAUGCCUGUGGGUCACAACAGAAAAUACAAUUUAUCUUCAGAGGACCUGUUUCUAUUCUGAUUGGGAUUGUAUGCGGUGCCUUUCAGACAAGGGAGGCUACGAUGGCAAUUAUUCGAGGGUAACCGCUUAUCACUAUUGUCGGCCAAAAGUUCACGAUGAGGCUACAGUGUUAUACAAACUUAAGGGUACCCCCAUACGAUUUGGUCCGGCACACAACCACAGAGGUGCCUUAGCUCUUCGUCCUGUUGAAACCACCCCCCUGAGGACAUCUCCUGAAUAUGUGGAACACCCUUACGACGAGGAGUCUAAUACAGGUCGGGGGGCAGCUCGGCAGGAAGGGCUCGAAGUUCAGCUCAUUGCAGGCCAGUACGACGCUUGGCAAUGUAUGGAUAACUGGACGGUUCACAUCAGUCGACUAUGUACCGGGGUACAGUCCCCACCGUUGGCUCGUAGAGUGUCCUAUUUAGGAAUCAGUCCGGACAUUUGGGGCCAUAGUUUUAGUUUCCGCUACACUUGGCCCUUCAAAAGCUCGACCGGGGCGAUUGAUUGGCGCCACGGGUUGUACUCCUUCCUGUCUACCUUAAUAAUCGGUUGCUGCUCACGAGCUCUGACUGAGAGAUCACGUCGAGGACCAGCCCGUUAA";
            string expected = "MPMGLVWHKQGPLERISIRGVIGVRSGYNETIRRNWVMLVSKSALFVSTCCQCNPPYLTCYKQLKSPDVTRFARAHDMDHFRDHTHMAGRTNLEQTFCAQPVHLTMDLEYYQDPPVAYVLQYMVPRRMPALSVITNPPNQELHSLWSYLHLSVSQISRGRLTCVTLLLVIWDNHLLRYSYLCQAFEMVMSSNIRSLGAETLKVLTDPNRECASVTIPHLTQQSALSIRAPPRYLHTELRRCRTSLARTTLRFVSSCTSRFGASIPSFHELYCSLYTGLSFVHISRSIKNAVDQTCGKPIDNTQIVISFQSSRLHRCSLIARLQADRTFFLYKRGWCRGVCSIVIYLTGICYVKLCCYATGISTVSDSYKPRGYVHLQVFLFRHVTVYVIRTYENSVTGASRYQADYIMGRDLPTLVDSLPVVSYSWNQYLLLVCGPRRGDLLWPLLGRDPTVKFKPRTYRQRRSFTAFMTRKCVSLWAQLTARRRTSRKLIRALTHGWKMRMSRHAVVYTGIEGILGSKTELTVPTHPLYILLPCSGPSCVLMTGVRPFWSSSQDSRWSTTIRPEGGVLSVSCRQIANSTFSVLESLGLFIEVRHGHGKIPLYRSSTCSNCSHVCQSNEWTAWFLNSPAAGPNQCQIVYNTKYCIAGYAPSPLLSTARAASYRYKCSWYFLLLFLTCSLISGQFLEPGIMKCALGMLKVGSECNLDEPHVSSTAGYHGHVNTLCCIIAVSPPLSLRVKRPAGRVVFEIGPACDQDGGQLASNIVITRNPTGWHPWTSRVSALTIRSLQCPQCVNEGWHWSHPEPGHMGRYLHPNSIPYEACVSTYTIAEKFSVHFSTRSIDWECPVRILATYTRQLELRCLTGTCFSYWNTPKIPTVTPHVARKSKETDRNASIQPPMAATVPATGLSKQWIVIQSSRDLRTDDEMDTTFQTHRYGFGISTTECYTLPSLYSAPKSRNTHTCHASPSNDSWSCFKLCGVSIRDRKRRVWTASRCNDLASEKGGRHGTPPTRIPPSLPRKTTQTGAANLAYWSQFINGDKSPDGYTKHVGTDAIYAARVQAHSLTCSRLSRIVCSGSAVGWDGDDFRESDVDRECEVTSGPRPQRFCIYIMSSVSTDGGVCRKIWIFHLLSNSCRSLFSPSSDREQKLGIMAQAEVFCHALSKFLFRLCSISSLRALSSYPRSGRVLSIVTNIPQWDSRPRLRIDYHQFTPIFHNPPGKLLQPNDTCLGCSNCLQSSPLLFTRAGLLGLFRTMAPITYGCMTRSAWSSLQRGGPLDVILSTSRATLFKRCVRSRFRNPCGINVLTNFHEPRALLFSAPRGNPSYTVMFASTAATCERLALARHRTYSRFLFERTRKHALKGRTSDQDGLSAGVTWKGVNSGRYLYSTCRHGLRRLIRTMMKQIGRLFYRTRCTLLVHRDPGAEHYLSTPIVKQGRGYHRRTRLPFLLFYDQCPARGVPLLVRLVTNALSPYITMSYACQTCNRAQESVSLLPSTCAYPYRSINSDPVGAWKVGSLNTHYICKLPKVSPPRWASSQHIASTVVRSDTHSRDESINGHGRLKVAIIVCRVTARSVAEVGGYWKKRHRSPRSIFDMKNRDTWHDGGRGPTSSSDHYYPITGIPVLSTPSRRPQARLTIKSSSIWAYRDQPMLTQPCITYPVLLLQWAVASRMSAISGIRLAGLWGAANAPVQLFLYLSQCAHLGTYCVYHEMGPMYMVLPNQRRLRRFSLGGQILCYPFRTTYSRAVVRPYPDTPHGYSGNGVNTGFLRRWPGPVPCLWVKGVAHHASKLRWKTLACTALSLPALPSHPGLLKSLLKRGEVKKLEAARTVLRALQPHRTNGGYTKELGQEEEHWFRYSFGGTIGALYDRKVAVWYRSSQVIYNGPIYRLVQEPSCASPITLFSRVLLNVVSYLEITVKRHLFTPNTIPGRGLSHRKGALFFLIPPYPSSFELIRKPYRRTQPAFPGSLGDDPLLPCICRSKLRLRRITVYGSSRPVEDVGTTESRCPKTFYSLPLALCWCRSTGVGIYGYLEPDYPYYTRRAPAIPIWRHPVYRAFVQLTLAQATVPGPKPFLLLPVGTGPFIKIFGRHGFAVLAAVRYCTLCSDPIRGPDSRLSLVTESARQVLTPVEYTDLQLCSLGSIQKSKRLSRLSLTKFDELAARTNARFVRPVGYQQVKPRVIETETPYRTYRPNPPTDVDEEHKHRLLAVKRYIGFRGTLSDRPQAALIKGECKASRSGHLLEIHRGSERLDCRITLTLLKDPLGYFYGVDSQPPAVHLSYDAFSQLIKNRCFRNYKSRAFQGTNVARKSIWDRSSDPGLREPLCWTLDKRWAKYTHQPQPCQPVPTSNVFECRRGQTEVRHNALVYSPTFANFRVSADSTVLLPATGSLIVPSGRSERNPDYHAPCSLCPNLNSGLPGRIVGTPRTDKRRAVTQLRVSEASCPQRNRSPDSCALGKCPLPLDNSTSRFCGPCAHYSDLTTFLRLLFSQTSPDGRFFHQIIAVLKVGKCIRATAALNPASYLASRAAQSGNIVEGHRGPILSRSMLVKRNWRCHGSLPCNRTTGETLTILGIMYGSVGGSRILQAVVHSVARRAPQATDGFGHLLTKPFVRARTSEEDKGWVVTSIFPRRNAYDPLGKVGPTWRTPGYRMIPQADKWGIGPTNGDTANCVGGCLRVICPSRASKHHAAICLSGTSFARVGKGREEKVRHTESQAWTRFNLGNARRGGVDQILRKLTIDRGGLRGHARISTKYRPRLQRNATMAGKDTRHYLSAFKFELATRAYRSALILWELIQRSTRTLAYSTGRIELPTTVPHSSPGVVLLGVVNLRGEATQHALFSYLDSTTPYECLDIAPNSWFTSESSWPRSDAILSYRLVNDPAELLSAGPSTLSHPLCGSIGRSLTRGDRLTKDYYVVHNLSSSLTRLDRHSSATRSLRPCLWVTTENTIYLQRTCFYSDWDCMRCLSDKGGYDGNYSRVTAYHYCRPKVHDEATVLYKLKGTPIRFGPAHNHRGALALRPVETTPLRTSPEYVEHPYDEESNTGRGAARQEGLEVQLIAGQYDAWQCMDNWTVHISRLCTGVQSPPLARRVSYLGISPDIWGHSFSFRYTWPFKSSTGAIDWRHGLYSFLSTLIIGCCSRALTERSRRGPAR";

            // act
            string actual = BioInf_PA2.ProteinTranslation(pattern);

            // assert
            Assert.AreEqual(expected, actual);
        }


        [TestMethod]
        public void PeptideEncoding_TestPattern1() {
            // arrange
            string DNA = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA";
            string peptide = "MA";
            string expected = "ATGGCC GGCCAT ATGGCC";

            // act
            string actual = BioInf_PA2.PeptideEncoding(DNA, peptide);

            // assert
            Assert.AreEqual(expected, actual);
        }

        [TestMethod]
        public void PeptideEncoding_TestPattern2() {
            // arrange
            string DNA = "TTCGAATCCGTATTGCTTTGGTTGGTTGAGAGAGCGACCTACATTGCTAGTCAGAATAGGTGATTCACACAAAGCTTAGCACCTGGGCAGCACCCCGTGATGTAAACCTATGGGAACTAAGGGAGTCCTGCGGTTTTAGCCAGCAAGCGAGCCGGCAGGAACACTCATACCATCGGACGCGTTTGACGCCTCCCCGGAAAGGAAGTATTTGAGCCTCATTATTACGTATTGCCCGTTAGTCGACAAATCAAGCCCTCGTACGCAGCTTATTCGTACGACGTGGAGGCGTTCCCACGGGCCTAACACGATTGGAACACCACCATAGTAGTGTGGTTCAAATACCTCCTTTGGAGATCTAGAGCTTCACTCTGATTCTAGAGGCAACTTTACAATCGCTCTACGAAATTGTATGGACATCATCAACCGGATATTCTGGGGCGGTAGAATTTCTTTTGTTCGAATCGCTCTAGGCCAGGATCAAATTAATTGAATTGCGGACTCAAGGATCGCGATAGCCGACACATCGGACGCTGTAGAAAGCCAGTCTCTGGATTTAATCCACCCTCTATGTTTGACAAAGCACTAAAACGGGATAGTTTCGGGTGGTATAAGTTTCCCAAGACGATTGCATCGCAATTCATCAACAACCATGAACTTACTGTTTTAGTACTTCCACACACCTTGTTAAATTACGCCTTTACTTCATGTTGCGGTGTGTGTTAGATAGTGTGCAGCTACAAGTCTACCGCCATCGCAGCTCGGGATACCGGCAGATGAGATGGTCCTGAGCTCGTACCGGACTCAAACTTTTTCCTTTACTACCTAGGAATCGCCCATGCGAATTTGTCGGACACACACCATTACATTAACGTCACAACAGCTACTGTTAGAATTTTGCTCTTGCAAATCCTGGAAAGAGTTAAAAAAACTCTTCCGCGCGCCAATAGGGTAAATAATAGATAGCCAGACGGCTGTAAGAGGTGATGACATTTGCAACAATCATGCTGTCGCATCTTCCGCAAGTTCATGTCGCGCCTAGGCAATGGATCTGCGAATGGGGGCCACGGGGTATGAACTACGGAATTCTAAGAAAGTTGCCATCCAGAGTTAAGGGTTTGAGGCTAGTTGCATCGCTGGTAACGAACTACCTCATTACTTGGACGCGCAGTGTGACTTCACTCCTGTATAGCGATGATGCCAAGCAGGAATTAGCAAATCTGAAGAGCGTTTCCAAACTGGCCACTTGGACTGACACCTATCGCGGGGGATTTCAGCGCGTGTCGCTCTCACATGAGAGCTGCCGTCAGGAGCGGTAGAGTTTAGAGAGGAATGCGACAAACTCCCTATTCACCTCTCTGGTGATGTAAGGATATTTACGCTTAGTTCTATGCCAGGCTTAGGGCCTCTCGGAACTTTGGTGAGTCCTTATTAATTGATGCTACCTCTCCCTTACCTTCGCCCCAAGTCACGTAGAAGTACTCAATCCTGCTACATGATAATCAAATATTTCCAACGTTGGGAAATCGGTGACATCACATACTAGTTAAGAAACCACTGTCAGTGAACTTATATCCGGGGGAGAAAATCTACTAACTTACATACGCTGTGCGAGCAGTTTTCATTATAAGAAAATATACTCCCGAGGTACCGCATCAAGCACGACATTCCCGGAGAGCATAACATTTCGGTGCACCTGCTTTTGTGCGCTTGCTTGCGGTTATTTATAAACTACGCACAAGGCGCAAACCGCAGTGCGCATGTTTTCTCCGCCTGGCTAGAACTCGACATTCTCGTCAACGCCAATCTATGTGAGAGGATTTAGACCTCTGTGAAAACGAGTCCCTCTATAGAATAAATACCCAGATGCCAATGGGGGTTCTATCCGATGGCAGTGCATGGAGTGGTGGCTCCAGATTAAGGATGAGGAGAGGTAAAGATAACAGTTCGGTCGCCACGACGCGTTGCCAATCGAAATATCAGTACTAAAAGGCCCACCGCTCCGCTTTAGTCCGACTTTACATCCTGTGGAAATTGTCGAACGGAGGCTACATCGGGCTATATGAGTGTGAAAACCTATACTTCTCGCGTCGTTACTCAGTGCCGGTCTCCTGTTTCCCCCAGTCTTACGTACCCTTATTGATATTTGCTTCACGTTGAAACGTCCTAACGCAGCGTAAAGAGGTGTTTGAACCTCATTACTATAAAATCGCGATCGAAGGTAGACTACATACGCAAACGCCGAAACCCTCAGTTGGCCTTGTTGCAAGTATGGAACGTTGTAAAATTTTTCCTAGACGTTGAGCTATCGGTACAAGGTCGTTAGCGTCCTTACCCTTCACTTATATGCCCGACAAAACGCGGGTCCTAGTGCAGTGGTGGGAGCTTGGAATCCCGCAATACAAGGACAACCTGTATCTCGTTCGGCGTTCCGCGATCACTCGATCCCGAACCACTCCAAGCCTGGTTGATCAGCAAAAGCGGAAGGATGGATAAAGGGCTACTGGTTAATGGATGTAAACTTCCAATGATGAAATCCTGGAAACGAGGGATCGGGTTACGGTGGCGAACGGGGTACGGCAACGTGGCTATCTAGAGCCCGACGTTACGACTCATGTACATGCTGCTACGTGGTTGAAGCTGACGTTCAGATGAAGCAGTACTGAGTCCTAGGGCTTACTACTACTCCAATAGGTCTGGCCGGCCAGATACAAAAGTTCGTGGCGGCTCACCCCCTTTCTGGCGGGTGTAGCTTGCTGACCGGTTTGCTCGATAACACAGGCTAGCGAATAGTAATGAGGTTCGAAAACCTCTTTCCAACGACTGAAAGGGTCTACACGAACTATCTACATTTCCCCGCCCATGTCCTTCCGTCTGGTTGCTTCTGGAGATCCTTTCGCATTATACCGCAGCGTAGTGGCTCTGGCATATATGAAAAAATCCTTCTGTGGGTATTTGTGCCATCACTTATTGTTCGTACCGATATGGGATTACAAGTGCGATGTGATAATAAGCGAAGAAGCCAACATGTTACACTGTTCATGCGCTCCGGGTAATGTGCGGGCACCATGCTCAGTTCCCGCTCGCAGTTGTCACTGTCCCTGTTTCGGCACCATAATCAACATTTCCACGGCCACGCTGGTGAATAACCGAGGATACCGAAGTACAGCAAGAATGAGAGCGGGACTCCTCCATCTGCTTGTAATACGCCTTCAAGATAGTCCATAAAACGGTCGGGGTCTTGTTTCGGACTAGCCGCTTTGAAACGGTGCATAGTTGTGTCAAGTGTGGACATTGGCTTTCTATCCTCGTCAGCGATCCTCGGAAAGACTCGGGCAGTCGCCCCGAATCGTAATTAGGTAGTAGTGCGGCTCAAAAACTTCCTTCGACCTAACCGCTATAATGTTCGTAGATATAAATTTCGTTTCAGTATTAACAGGCGCACCGTATATATACGGAATGGTGTCGCCCCATTAGCTGCTCGCCAATATTTATCTAAGACCGCGCGCGTCTAGCGCCTTTAGTAGTTGCACCCGAGTATAGTAATGGGGTTCGAAGACTTCCTTCGCAAGGCTGCCATACTGTATCACAAGTACTGACGGAGCCCCGGAGGAGTGCAGGATACGGCAAAGGAGACCATTACCGGGGCATGAGTCCAAGTTAGCCCGTTAGGTGAAGGACGCTGATACAATAGTGAATCCGTTACTGAAAGGTTTAGAAGACCGGGGGGCTCGCACTAGGTCCAAATATTATGAACCCTACTCCTGCAACTGAATTGGCCGTCCAGGCGATATTTAAAAGGGGTTACTAGCAGGTTCATCGGAGCCCGTACTCCTTCCGGGCATAGTCGTTCGACGGGTAGAAATTCATCCAGTCGTGCCGGATACCCCGAGAATACCCCTATTTTTTGATCCTTCACCATCATCGTCCGCGGACTCATCTAAGTACCTCAGACCGAAACTGTTATCGTAGCGAAGAGCGAACTCGAATGACATCGCTTGTCCAACAGGGAAAATATGTAAAGTATATGCAGATTATTATAGGAAGATCACAAACTCCATCGCGCCTAGGCCAAAGACTTGCCAGAACAACATCTCTTCCAGAGCAAGGAAGTGTTTGAACCTCACTATTATCGAGAGAAGTCCCATGAATTTATAATAGTGAGGCTCAAAAACTTCCTTCATCGTCGGGCGCTGGGGCGAGCTAGGCTTCCCTAGCCGTCATTACTGTACCCACGCCAAATATTATCGAGTATGACTACGAAGCCTTCACAAGGGGCTGCAGTAACAGACTAACTCACGCACACCGCAACTACTATCCTAAATTGAGTAAGAAAACCCTCGACTATAGCCCTTGAATTAAATTCGCTATTTAAGGAAGACCGCGCTTCCGCTCGCCCGCGTATAGTTTACCAGTCCCAAGAAACATGTGTGGCCAGCCTACCTGAAGAGCTGTGTAAAGATAGATTCTGACATCCTCAAAAAGAAGTTTTCGAGCCGCACTACTACGCACGGAGCTCCGTTATTCAAGGCATGTCGAAGTACAGCGTGGTGGCCTCTCCTGTTCTCCACCCCAGCTAAACCCACCTTGTTCGAATTCGCGCAACTGTATATGACATGAACACTTACAGGGGTTAGAAGTTACTGCAACCAAGATAGAGCTCGTCGAAGTAATAGTGCGGTTCAAAAACTTCCTTCAATTGGTCTCATCACTTAAATTTAAGAGCTATTGTGAGTACAGGTACGGATGCGGCTTCAGTGGATCTTCAGCATTCATTCCTTGTAGTAATGGGGTTCGAAGACTTCCTTGCCAGGGTACCAAACAAGTCTTGCGCATCCTCCTCCCTAAGGAGGTATTTGAACCCCACTATTACCCACGATAGAACATGCAGGGTTTGATAGTGGAACACCTTTTAGAATCTGGGGATAAATTCCCAGGACTAATGTATGGCTGTAGTAATGAGGTTCAAAAACCTCCTTTTCAGGTGGATCGCAGGCCGTGCTGCCTCACAAGCTGGGACGCCGTCCACGGTATAGCCGGCGTCGGCAGTTACTGTGAAATAGCGGAAACTCGATCCCAATATATCATCTTACGTTTGGCGCCCAATAGTCGCCCAGTACCCGTTGACAGTTCTTTAACTCGGCTTAGAACTACTAGACAGGTTCAACCGAACCTTGCCCTAGTTCCCACTCCCGTAATTCATTTGGGTTTGCATCTAGAATACTGGAGGGTGCATAGACGAAACGTGTACGTCGGAGAAAACGTAAGAAAATGACCTAGACTCATAGTAGTGAGGTTCGAAGACTTCCTTTCAGTGAAATCGATCCACCACTCGCCGCGAAGAGATAATAGCATAGAGCACAAGTGCGCGAGTAGAGAAAAAGGCTATCCCAACCGGGCACGTCCTTCGTGTTTGGCGTTTACATACGGCACCCCGTTTCTGCACGTTAACCGTCTAGTATCCAACGGTGGATGGCGGACGCTAGACTATAGATATGAGATATCGAGACCTGGAGCTGGGTGTGGCTGCAGCCCGGGTCATTGCGGGCTGTGAATTCAAGGGCATGTAAACAAGCGTATATCGAACAGTGGATGGGCACCTGCAATACTCACGGTAGAGTTAGCTCACAGGATTCACGTTGAGGACTATGAGTCCCTCTTCGCTAGCAGTCTGGGGGGATATGGAGTTTAATAGCTTGACGTAGTAATGGGGCTCAAACACCTCTTTGTGTGAGCACAGCTACTTGCATTAAGAGATTCTAAACAGCGATCATCTCGGCTATTTCGGGCCAGCCTTTTCGGCAGGATGTTATGTAGCATTTCTGGAAGCTTCCCCCTCGAATCTACTAGTGGTGAGAAGATGCCCCACCGATATTACTCTTTAATCTTGAGAAACCTAAAACCGATCTGACCTCAGACGGGCGGCTCCCACCCGAGGATAAACTCGTCAATAATAATGTGGCTCGAACACTTCTTTTCTCACTAGGCTTTTACGACACGCCACATGTATTTAAGCATCTACCTAACTTGTGTCTGCTGATATACAGCGCATTCTACGCCCAACCTACCAATTACTTCAACGTAGTGCGTGGCTAAAATTCAGGGGAGCTTCATCTCTGTCTTAATTTGAAGGTTCTTCCGGGGCGTTTGGGAATCTTCGTGCCTTTTGCGAGGTTAAGGTATCAAAGAAGTTTTCGAACCACATTATTACCGCCTTAAGCACCGGCGCATCCTGCTCGTGACAACTCTACCCTGCCCTGATAAAGGCACTGAACGTTCCAGAGAGTGCATCATTGACACGCGAGCAGGCCACAGTAGCCACAAGACGTATGGGTGATTATAGAATTGGTGGAGGTGTTGTTAACGATCAGGAGGACATTAGTGGGAGTTAGGAAAGACCCTATGTTCTCTCTATCGCGGACTTGTAACTTGACAAGCAAAAGGGTAAGAGAGCTGCACACCGAAGCAGGCCCTTCCTATACCTGTTTTTCCTACGCGTAGAGAGGAATCCAGAAAGGTGATAATTGGCATTCGATGAAAAAACAGTGTGCCACTGACTTAGTTCTATATGTGAAGAGCCTGTTAGCACGTGACGGCGGCCTTGGTATAGAGCCCTTAATGGTCTCCATCGCGTAGTAATGGGGCTCGAAAACCTCCTTACTTGGGATTGCGTGGCCTCCTTGTGAGTCATACACAAGGCTTAGGGCTATGGGGCGATACACTCCTTTTCGCGGCGCATGGGGCGGTGATGCCTACATAGTAGTAGTGACTGCCTTTCTGGGGGGCTATTTGTGGATGACCAACACCTGACCAGCGATGCAATCGCTAGGGGAGGTACACCTCTCATATGTTACAACAATCACCGAATTGTGTTTCGAATTCGAATCAAGTTTGCGGTGTCGACCAGATCTGGTCTTGCTGCCATACCGGGTTCGCCGCCTCCGGTGGATAGAACTGCATCTTAAGACATCTGGACCCAGCGGTAAGTAGCGGGAAGAGTTTAGAGTCATTCGTACAACTACAGGCTAAGGGCTTACTGGGGAGTTGTTGTAGGGCATAAAGATCGCCCCATGACTTTTCGTACTTTCCCCGATAGTTCACTCGCAGCGAGCTGCGGCTGGGCTTCGCCACACGAGTACGGGCAACATTTATCTCCTCTAATCACTGGGCACCGCGCGAGGAAATAGAAAACCCTAATCAGTGCTCATGGGCGCATCTATTGGTCTCCGCATGCACGATGCCGCGGAGTGCTTAGTTGTCCCTGCATAATCTTCGTAGATGTATAAGAGATTACCTATTTATTCGGTTTCGGTTCTAGACGTACCTTGCCGCATGAGTATAGGCTAATGAACTGAGTTGGCGCCAGAGGGAAAGGCATAATAATGCGGCTCGAATACTTCCTTAAGGAAGTATTCGAACCACATTACTAT";
            string peptide = "KEVFEPHYY";
            string expected = "AAGGAAGTATTTGAGCCTCATTATTAC AAAGAGGTGTTTGAACCTCATTACTAT AAGGAGGTATTTGAACCCCACTATTAC AAAGAAGTTTTCGAACCACATTATTAC AAGGAAGTGTTTGAACCTCACTATTAT AAAGAAGTTTTCGAGCCGCACTACTAC AAGGAAGTATTCGAACCACATTACTAT ATAATAATGCGGCTCGAATACTTCCTT GTAGTAATGGGGCTCGAAAACCTCCTT GTAGTAATGAGGTTCAAAAACCTCCTT GTAGTAATGGGGTTCGAAGACTTCCTT ATAATAGTGAGGCTCAAAAACTTCCTT ATAGTAATGGGGTTCGAAGACTTCCTT GTAGTAGTGCGGCTCAAAAACTTCCTT ATAGTAATGAGGTTCGAAAACCTCTTT ATAATAATGTGGCTCGAACACTTCTTT GTAGTAATGGGGCTCAAACACCTCTTT ATAGTAGTGAGGTTCGAAGACTTCCTT GTAATAGTGCGGTTCAAAAACTTCCTT ATAGTAGTGTGGTTCAAATACCTCCTT";

            // act
            string actual = BioInf_PA2.PeptideEncoding(DNA, peptide);

            // assert
            Assert.AreEqual(expected, actual);
        }
    }
}