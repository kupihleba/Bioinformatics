package main

import (
	. "Bioinformatics/lab1/needleman-wunsch"
	"Bioinformatics/lab1/utils"
	"testing"
)

type test struct {
	seq1 string
	seq2 string
	res1 string
	res2 string
}

var testsDefault = []test{
	{"A", "A", "A", "A"},
	{"TAGA", "TCGA", "TAGA", "TCGA"},
	{"ATA", "AA", "ATA", "A-A"},
	{"TAGTT", "TATT", "TAGTT", "TA-TT"},
	{"GCATGCV", "GATTACA", "GCA-TGCV", "G-ATTACA"},
	{"AACTTATAGGGTAG", "CTTFATAGVGTA", "AACTT-ATAGGGTAG", "--CTTFATAGVGTA-"},
}

func TestAll(t *testing.T) {
	engine := NewAlignEngine(
		ScoreDefault,
		-1,
	)
	for i, test := range testsDefault {
		res1, res2, _, err := engine.AlignSequences(test.seq1, test.seq2)
		checkTest(err, t)
		if test.res1 != res1 || test.res2 != res2 {
			t.Error(
				"\nFor sequences:\n",
				"\""+test.seq1+"\"\n",
				"\""+test.seq2+"\"\n",
				"Expected:\n",
				"\""+test.res1+"\"\n",
				"\""+test.res2+"\"\n",
				"Got:\n",
				"\""+res1+"\"\n",
				"\""+res2+"\"\n",
			)
		} else {
			t.Log("TEST", i, "OK")
		}
	}
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
