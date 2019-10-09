package needleman_wunsch

import (
	"Bioinformatics/lab1/utils"
	"fmt"
	"io/ioutil"
	"log"
	"strings"
)



type AlignEngine struct {
	ScoreFunc ScoreFuncType
	ScoreGap  int
	GapChar byte
}

func NewAlignEngine(scoreFunc ScoreFuncType, gapPenalty int) AlignEngine {
	return AlignEngine{
		ScoreFunc:scoreFunc,
		ScoreGap:gapPenalty,
		GapChar: '-',
	}
}

func init() {
	log.SetOutput(ioutil.Discard)
}

// Bioinformatics algorithm
// to align protein or nucleotide sequences
func (engine *AlignEngine) AlignSequences(seq1 string, seq2 string) (string, string, int, error) {
	table := make([][]int, len(seq2)+1)
	for i := range table {
		table[i] = make([]int, len(seq1)+1)
	}

	// Init first row and column
	for i := range table[0] {
		table[0][i] = i * engine.ScoreGap
	}
	for j, row := range table {
		row[0] = j * engine.ScoreGap
	}

	// Fill the table
	for i := 1; i < len(table[0]); i++ {
		for j := 1; j < len(table); j++ {
			score, err := engine.ScoreFunc(seq1[i-1], seq2[j-1])
			if err != nil {
				return "", "", 0, err
			}
			match := table[j-1][i-1] + score
			del := table[j][i-1] + engine.ScoreGap
			insert := table[j-1][i] + engine.ScoreGap
			table[j][i] = utils.Max(match, del, insert)
		}
	}
	printMatrix(table)

	return engine.findAlign(seq1, seq2, table)
}

// Find align for both sequences with the given weight table
func (engine *AlignEngine) findAlign(seq1 string, seq2 string, table [][]int) (string, string, int, error) {
	i := len(table[0]) - 1
	j := len(table) - 1
	resScore := 0
	var sbSeq1, sbSeq2 strings.Builder

	for i > 0 || j > 0 {
		if i > 0 && j > 0 {
			score, err := engine.ScoreFunc(seq1[i-1], seq2[j-1])
			if err != nil {
				return "", "", 0, err
			}
			if table[j][i] == table[j-1][i-1]+score {
				sbSeq1.WriteByte(seq1[i-1])
				sbSeq2.WriteByte(seq2[j-1])
				i--
				j--
				resScore += score
				continue
			}
		}
		if i > 0 && table[j][i] == table[j][i-1]+engine.ScoreGap {
			sbSeq1.WriteByte(seq1[i-1])
			sbSeq2.WriteByte(engine.GapChar)
			i--
		} else {
			sbSeq1.WriteByte(engine.GapChar)
			sbSeq2.WriteByte(seq2[j-1])
			j--
		}
		resScore += engine.ScoreGap
	}
	return utils.Reverse(sbSeq1.String()), utils.Reverse(sbSeq2.String()), resScore, nil
}

// Only for debug purposes
func printMatrix(matrix [][]int) {
	var sb strings.Builder
	for _, row := range matrix {
		//_, err := fmt.Fprintf(&sb, "%v\n", row)
		for _, i := range row {
			sb.WriteString(fmt.Sprintf("%3d,\t", i))
		}
		sb.WriteString("\n")
	}
	_ = log.Output(0, sb.String())
}
