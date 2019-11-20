package main

import (
	. "Bioinformatics/Sequence_alignment/algoritm"
	"Bioinformatics/Sequence_alignment/utils"
	"fmt"
	"testing"
)

type test struct {
	seq1 string
	seq2 string
	res1 string
	res2 string
}

var testsNeedlemanWunsch = []test{
	{"A", "A", "A", "A"},
	{"TAGA", "TCGA", "TAGA", "TCGA"},
	{"ATA", "AA", "ATA", "A-A"},
	{"TAGTT", "TATT", "TAGTT", "TA-TT"},
	{"GCATGCV", "GATTACA", "GCA-TGCV", "G-ATTACA"},
	{"AACTTATAGGGTAG", "CTTFATAGVGTA", "AACTT-ATAGGGTAG", "--CTTFATAGVGTA-"},
}

var testsSmithWaterman = []test{
	{"ATA", "ATAT", "ATA", "ATA"},
	{"CTCTGAG", "TGTCAGT", "TCTGAG", "TC--AG"},
	{"CTCTGAGG", "TGTCAGTA", "TCTGAG", "TC--AG"},
}

func TestNeedlemanWunsch(t *testing.T) {
	engine := NewAlignEngine(
		ScoreDefault,
		-1,
	)
	for i, test := range testsNeedlemanWunsch {
		res1, res2, _ := engine.NeedlemanWunsch(test.seq1, test.seq2)
		//checkTest(err, t)
		if test.res1 != res1 || test.res2 != res2 {
			t.Error(format(test, res1, res2))
		} else {
			t.Log("TEST", i, "OK")
		}
	}
}

func TestGapDynamics(t *testing.T) {
	engine := NewAlignEngineDyn(
		ScoreDefault,
		func(gapsInRow int) int {
			return -2 - gapsInRow*2
		})
	test := test{
		"AAAA", "AAAAAAAAAAAA",
		"-A-A-A-A----", "AAAAAAAAAAAA",
	}
	res1, res2, _ := engine.NeedlemanWunsch("AAAA", "AAAAAAAAAAAA")
	exp1, exp2 := "-A-A-A-A----", "AAAAAAAAAAAA"
	if res1 != exp1 || res2 != exp2 {
		t.Error(format(test, res1, res2))
	} else {
		t.Log("OK")
	}
}

func TestHirschberg(t *testing.T) {
	engine := NewAlignEngine(
		func(a byte, b byte) (i int, e error) {
			if a == b {
				return +2, nil
			} else {
				return -1, nil
			}
		},
		-2,
	)

	res1A, res2A, _ := engine.NeedlemanWunsch("AGTACGCA", "TATGC")
	//checkTest(err, t)
	//res1B, res2B, err := engine.Hirschberg( "G", "GC")

	res1B, res2B := engine.Hirschberg("AGTACGCA", "TATGC")
	//checkTest(err, t)

	if res1A != res1B || res2A != res2B {
		t.Error("Needleman-Wunsch and Hirschberg results mismatch!")
		t.Errorf("\n%v\t%v\n%v\t%v", res1A, res2A, res1B, res2B)
	} else {
		t.Log("OK")
	}
}

func TestSmithWaterman(t *testing.T) {
	engine := NewAlignEngine(
		func(a byte, b byte) (i int, e error) {
			if a == b {
				return +2, nil
			} else {
				return -2, nil
			}
		},
		-1,
	)
	for i, test := range testsSmithWaterman {
		res1, res2, _ := engine.SmithWaterman(test.seq1, test.seq2)
		//checkTest(err, t)
		if test.res1 != res1 || test.res2 != res2 {
			t.Error(format(test, res1, res2))
		} else {
			t.Log("TEST", i, "OK")
		}
	}
}

func format(test test, seq_res1, seq_res2 string) string {
	return fmt.Sprint(
		"\nFor sequences:\n",
		"\""+test.seq1+"\"\n",
		"\""+test.seq2+"\"\n",
		"Expected:\n",
		"\""+test.res1+"\"\n",
		"\""+test.res2+"\"\n",
		"Got:\n",
		"\""+seq_res1+"\"\n",
		"\""+seq_res2+"\"\n",
	)
}

func checkTest(err error, t *testing.T) {
	if err != nil {
		t.Error(err)
	}
}

//func TestScore(t *testing.T) {
//	engine := NewAlignEngine(
//		func(a byte, b byte) (i int, e error) {
//			if a == b {
//				return +2, nil
//			} else {
//				return -1, nil
//			}
//		},
//		-2,
//	)
//	seq1 := "AGTACGCA"
//	seq2 := "TATGC"
//	//arr, err := engine.NWScore( "TATGC", "AGTACGCA")
//	//checkTest(err, t)
//	arr, err := engine.calcGridScorePart(false, seq1, seq2)
//	seq1 = "CGCA"
//	correct :=  []int {-16, -12, -8, -7, -3, 1}
//	for i:= 0; i < len(correct); i++ {
//		if arr[i] != correct[i] {
//			t.Error("calcGridScorePart FAILED")
//		}
//	}
//	arr, err = engine.calcGridScorePart(true, seq1, seq2)
//	correctRev := []int {-3, -1,  1,  0, -4, -8}
//	for i:= 0; i < len(correct); i++ {
//		if arr[i] != correctRev[i] {
//			t.Error("calcGridScorePart FAILED")
//		}
//	}
//	checkTest(err, t)
//	fmt.Printf("%v\n", arr)
//}

func TestPrettify(t *testing.T) {
	prettyStr := utils.Prettify("AAAAAAAAAA", 3)
	if prettyStr != "AAA\nAAA\nAAA\nA" {
		t.Error("Prettify fails")
	}
}

func TestReverse(t *testing.T) {
	str := "abcd"
	if utils.ReverseStr(str) != "dcba" {
		t.Error("That is not a reverse")
	} else if str != "abcd" {
		t.Error("Destructive method!")
	} else {
		t.Log("OK")
	}
}
