// Bioinformatics algorithms
// to align protein or nucleotide sequences
package algoritm

import (
	"Bioinformatics/Sequence_alignment/utils"
	"fmt"
	"io/ioutil"
	"log"
	"strings"
)

// Matrix structure
//    Seq 1 i ->
//     _______
// S  |[, , ... ],
// e  |[, , ... ],
// q  |[, , ... ],
// 0  |[, , ... ],
//    |[, , ... ],
// j

type AlignEngine struct {
	ScoreFunc ScoreFuncType
	ScoreGap  ScoreGapType
	GapChar   byte
}

func init() {
	log.SetOutput(ioutil.Discard) // Comment out this line to view debug logs
}

func NewAlignEngine(scoreFunc ScoreFuncType, gapPenalty int) AlignEngine {
	return AlignEngine{
		ScoreFunc: scoreFunc,
		ScoreGap:  func(_ int) int { return gapPenalty },
		GapChar:   '-',
	}
}

// Memory optimised Needleman-Wunsch algorithm using Hirschberg trick
func (engine *AlignEngine) Hirschberg(seq1 string, seq2 string) (string, string, error) {
	// https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm Some ideas
	var resSeq1, resSeq2 strings.Builder
	if len(seq1) == 0 {
		for i := 0; i < len(seq2); i++ {
			resSeq1.WriteByte(engine.GapChar)
			resSeq2.WriteByte(seq2[i])
		}
	} else if len(seq2) == 0 {
		for i := 0; i < len(seq1); i++ {
			resSeq1.WriteByte(seq1[i])
			resSeq2.WriteByte(engine.GapChar)
		}
	} else if len(seq1) == 1 || len(seq2) == 1 {
		res1, res2, _, err := engine.NeedlemanWunsch(seq1, seq2)
		return res1, res2, err
	} else {
		mid := len(seq1) / 2

		leftScore, err := engine.calcGridScorePart(
			false,
			seq1[:mid],
			seq2,
		)
		if err != nil {
			return "", "", err
		}
		rightScore, err := engine.calcGridScorePart(
			true,
			seq1[mid:],
			seq2,
		)
		index := utils.SumAndMax(leftScore, rightScore)
		log.Printf("LEFT: %s %s\n", seq1[:mid], seq2[:index])
		leftRes1, leftRes2, err := engine.Hirschberg(seq1[:mid], seq2[:index])
		if err != nil {
			return "", "", err
		}
		log.Printf("RIGHT: %s %s\n", seq1[mid:], seq2[index:])
		rightRes1, rightRes2, err := engine.Hirschberg(seq1[mid:], seq2[index:])
		if err != nil {
			return "", "", err
		}
		return leftRes1 + rightRes1, leftRes2 + rightRes2, nil
	}

	return resSeq1.String(), resSeq2.String(), nil
}

func (engine *AlignEngine) calcGridScorePart(reverse bool, seq1 string, seq2 string) ([]int, error) {
	height := len(seq2) + 1
	width := len(seq1) + 1
	if reverse {
		seq1 = utils.ReverseStr(seq1)
		seq2 = utils.ReverseStr(seq2)
	}
	column := make([]int, height)
	gapInRow := 0
	column[0] = engine.ScoreGap(gapInRow)
	for j := 1; j < height; j++ {
		column[j] = column[j-1] + engine.ScoreGap(gapInRow)
		gapInRow++
	}
	log.Printf("%v\n", column)

	gapInRow1, gapInRow2 := 0, 0
	lastOverwrite := 0
	for i := 1; i < width; i++ {
		for j := 0; j < height; j++ {
			if j-1 >= 0 {
				score, err := engine.ScoreFunc(seq1[i-1], seq2[j-1])
				if err != nil {
					return nil, err
				}
				match := score + lastOverwrite
				del := column[j] + engine.ScoreGap(gapInRow1)
				insert := column[j-1] + engine.ScoreGap(gapInRow2)
				lastOverwrite = column[j]
				var index int
				column[j], index = utils.Max(match, del, insert)
				if index == 1 {
					gapInRow1++
					gapInRow2 = 0
				} else if index == 2 {
					gapInRow2++
					gapInRow1 = 0
				} else {
					gapInRow1, gapInRow2 = 0, 0
				}
			} else {
				lastOverwrite = column[j]
				if j == 0 { // TODO: FIX THIS
					column[j] = 0
					gapInRow = 0
				} else {
					column[j] = column[j-1] + engine.ScoreGap(gapInRow)
					gapInRow++
				}
			}
		}
		log.Printf("%d: %v\n", i, column)
	}
	if reverse {
		utils.ReverseArray(column)
	}
	return column, nil
}

func (engine *AlignEngine) NeedlemanWunsch(seq1 string, seq2 string) (string, string, int, error) {
	return engine.AlignSequences(seq1, seq2, false)
}

func (engine *AlignEngine) SmithWaterman(seq1 string, seq2 string) (string, string, int, error) {
	return engine.AlignSequences(seq1, seq2, true)
}

func (engine *AlignEngine) AlignSequences(seq1 string, seq2 string, local bool) (string, string, int, error) {
	table := make([][]int, len(seq2)+1)
	for i := range table {
		table[i] = make([]int, len(seq1)+1)
	}
	// We track max element of the matrix for local alignment
	iMax := 0
	jMax := 0

	// Init first row and column
	gapInRow := 0
	table[0][0] = 0
	for i := 1; i < len(table[0]); i++ {
		if local {
			table[0][i] = 0
		} else {
			table[0][i] = table[0][i-1] + engine.ScoreGap(gapInRow)
			gapInRow++
		}
	}
	gapInRow = 0
	table[0][0] = 0
	for j := 1; j < len(table); j++ {
		if local {
			table[j][0] = 0
		} else {
			table[j][0] = table[j-1][0] + engine.ScoreGap(gapInRow)
			gapInRow++
		}
	}

	// Fill the table
	gapInRow1, gapInRow2 := 0, 0
	for i := 1; i < len(table[0]); i++ {
		for j := 1; j < len(table); j++ {
			score, err := engine.ScoreFunc(seq1[i-1], seq2[j-1])
			if err != nil {
				return "", "", 0, err
			}
			match := table[j-1][i-1] + score
			del := table[j][i-1] + engine.ScoreGap(gapInRow1)
			insert := table[j-1][i] + engine.ScoreGap(gapInRow2)
			if local {
				val, index := utils.Max(match, del, insert, 0)
				if index == 1 {
					gapInRow1++
					gapInRow2 = 0
				} else if index == 2 {
					gapInRow2++
					gapInRow1 = 0
				} else {
					gapInRow1, gapInRow2 = 0, 0
				}
				if val > table[jMax][iMax] {
					jMax = j
					iMax = i
				}
				table[j][i] = val
			} else {
				var index int
				table[j][i], index = utils.Max(match, del, insert)
				if index == 1 {
					gapInRow1++
					gapInRow2 = 0
				} else if index == 2 {
					gapInRow2++
					gapInRow1 = 0
				} else {
					gapInRow1, gapInRow2 = 0, 0
				}
			}
		}
	}
	printMatrix(table)

	return engine.findAlign(seq1, seq2, table, local, iMax, jMax)
}

// Find align for both sequences with the given weight table
func (engine *AlignEngine) findAlign(seq1 string, seq2 string,
	table [][]int, local bool, iMax, jMax int) (string, string, int, error) {
	var i, j int
	if local {
		i = iMax
		j = jMax
	} else {
		i = len(table[0]) - 1
		j = len(table) - 1
	}
	var sbSeq1, sbSeq2 strings.Builder

	gapInRow1, gapInRow2 := 0, 0
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
				gapInRow1, gapInRow2 = 0, 0
				if local && table[j][i] == 0 {
					break
				}
				continue
			}
		}

		if i > 0 && (table[j][i] == table[j][i-1]+engine.ScoreGap(gapInRow2) || j == 0) {
			sbSeq1.WriteByte(seq1[i-1])
			sbSeq2.WriteByte(engine.GapChar)
			gapInRow2++
			gapInRow1 = 0
			i--
			if table[j][i] == 0 {
				break
			}
		} else {
			gapInRow1++
			gapInRow2 = 0
			sbSeq1.WriteByte(engine.GapChar)
			sbSeq2.WriteByte(seq2[j-1])
			j--
			if table[j][i] == 0 {
				break
			}
		}
	}
	resScore := table[len(table)-1][len(table[0])-1]
	return utils.ReverseStr(sbSeq1.String()), utils.ReverseStr(sbSeq2.String()), resScore, nil
}

// Only for debug purposes
func printMatrix(matrix [][]int) {
	var sb strings.Builder
	sb.WriteByte('\n')
	for _, row := range matrix {
		//_, err := fmt.Fprintf(&sb, "%v\n", row)
		for _, i := range row {
			sb.WriteString(fmt.Sprintf("%3d,\t", i))
		}
		sb.WriteByte('\n')
	}
	_ = log.Output(0, sb.String())
}
