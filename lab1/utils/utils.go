// Helpful methods that are not in std library
package utils

import "strings"

func Reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

func Max(a int, b int, c int) int {
	if a >= b && a >= c {
		return a
	} else if b >= a && b >= c {
		return b
	} else {
		return c
	}
}

// Inserts \n to fit lines into provided length
func Prettify(sequence string, length int) string {
	var sb strings.Builder
	for len(sequence) > length {
		sb.WriteString(sequence[:length])
		sb.WriteByte('\n')
		sequence = sequence[length:]
	}
	sb.WriteString(sequence)
	return sb.String()
}
