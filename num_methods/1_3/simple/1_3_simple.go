package main

import (
	"fmt"
	"math"
)

const (
	maxIterCount = 1000
	epsilon      = 1e-6
)

func SimpleIterations(A [][]float64, b []float64, epsilon float64) ([]float64, int, error) {
	n := len(A)

	for col := 0; col < n; col++ {
		pivotRow := col
		maxVal := math.Abs(A[col][col])
		for row := col + 1; row < n; row++ {
			if math.Abs(A[row][col]) > maxVal {
				maxVal = math.Abs(A[row][col])
				pivotRow = row
			}
		}

		if maxVal < 1e-12 {
			return nil, 0, fmt.Errorf("det == 0")
		}

		if pivotRow != col {
			A[col], A[pivotRow] = A[pivotRow], A[col]
			b[col], b[pivotRow] = b[pivotRow], b[col]
		}
	}

	// x = alpha*x + beta
	alpha := make([][]float64, n)
	beta := make([]float64, n)
	for row := 0; row < n; row++ {
		alpha[row] = make([]float64, n)

		beta[row] = b[row] / A[row][row]
		for col := 0; col < n; col++ {
			if row != col {
				alpha[row][col] = -A[row][col] / A[row][row]
			}
		}
	}

	// c
	normAlpha := 0.0
	for row := 0; row < n; row++ {
		rowSum := 0.0
		for col := 0; col < n; col++ {
			rowSum += math.Abs(alpha[row][col])
		}
		if rowSum > normAlpha {
			normAlpha = rowSum
		}
	}

	fmt.Printf("||alpha||_c = %.4f\n", normAlpha)

	// ||alpha|| / (1 - ||alpha||)
	coef := 1.0
	if normAlpha < 1 {
		coef = normAlpha / (1 - normAlpha)
	}

	xPrev := copyVector(beta)
	xNext := make([]float64, n)

	for k := 1; k <= maxIterCount; k++ {
		// x(k) = beta + alpha * x(k-1)
		for row := 0; row < n; row++ {
			sum := beta[row]
			for col := 0; col < n; col++ {
				sum += alpha[row][col] * xPrev[col]
			}
			xNext[row] = sum
		}

		// ||x(k) - x(k-1)||_c
		diffNorm := 0.0
		for row := 0; row < n; row++ {
			diff := math.Abs(xNext[row] - xPrev[row])
			if diff > diffNorm {
				diffNorm = diff
			}
		}

		epsK := coef * diffNorm

		if epsK <= epsilon {
			return xNext, k, nil
		}

		copy(xPrev, xNext)
	}

	return xPrev, maxIterCount, fmt.Errorf("max iterations reached")
}

func copyVector(b []float64) []float64 {
	copyB := make([]float64, len(b))
	copy(copyB, b)
	return copyB
}

func main() {
	A := [][]float64{
		{-23, -7, 5, 2},
		{-7, -21, 4, 9},
		{9, 5, -31, -8},
		{0, 1, -2, 10},
	}
	b := []float64{-26, -55, -58, -24}

	x, iters, err := SimpleIterations(A, b, epsilon)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
	}

	fmt.Printf("%d iterations:\n", iters)
	for row, val := range x {
		fmt.Printf("x_%d = %.4f\n", row+1, val)
	}
}
