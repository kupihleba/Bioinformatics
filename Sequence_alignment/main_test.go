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
		res1, res2, _, err := engine.NeedlemanWunsch(test.seq1, test.seq2)
		checkTest(err, t)
		if test.res1 != res1 || test.res2 != res2 {
			t.Error(format(test, res1, res2))
		} else {
			t.Log("TEST", i, "OK")
		}
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
		res1, res2, _, err := engine.SmithWaterman(test.seq1, test.seq2)
		checkTest(err, t)
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

func TestPrettify(t *testing.T) {
	prettyStr := utils.Prettify("AAAAAAAAAA", 3)
	if prettyStr != "AAA\nAAA\nAAA\nA" {
		t.Error("Prettify fails")
	}
}
