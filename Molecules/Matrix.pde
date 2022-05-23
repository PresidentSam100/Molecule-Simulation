public class Matrix {
    float[][] data = null;
    int rows = 0, cols = 0;

    Matrix(int rows, int cols) {
        data = new float[rows][cols];
        this.rows = rows;
        this.cols = cols;
    }

    Matrix(float[][] data) {
        this.data = data.clone();
        rows = this.data.length;
        cols = this.data[0].length;
    }

    boolean isSquare() {
        return rows == cols;
    }

    void display() {
        System.out.print("[");
        for (int row = 0; row < rows; ++row) {
            if (row != 0) {
                System.out.print(" ");
            }

            System.out.print("[");

            for (int col = 0; col < cols; ++col) {
                System.out.printf("%8.3f", data[row][col]);

                if (col != cols - 1) {
                    System.out.print(" ");
                }
            }

            System.out.print("]");

            if (row == rows - 1) {
                System.out.print("]");
            }

            System.out.println();
        }
    }

    Matrix transpose() {
        Matrix result = new Matrix(cols, rows);

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                result.data[col][row] = data[row][col];
            }
        }

        return result;
    }
    PVector vecConvertOut(){
      if (rows == 1 && cols == 3) {
        return new PVector(this.data[0][0],this.data[1][0],this.data[2][0]);}
      else if (rows ==  3 && cols == 1) {
        return new PVector(this.data[0][0],this.data[0][1],this.data[0][2]);} 
      println("ERROR CANNOT TRANSPOSE THIS MATRIX");
      return new PVector();
    }
    Matrix vecConvertIn(PVector in) {
      Matrix result = new Matrix(1,3);
      result.data[0][0] = in.x;
      result.data[0][1] = in.y;
      result.data[0][2] = in.z;
      return result;
    }
    Matrix multiply(Matrix matrix) {
      Matrix result = new Matrix(cols,rows);
      if (this.cols != matrix.rows) {
        println("ERROR: SIZE MISMATCH");
        return result;}
      for (int c = 0; c < matrix.cols; c++) {
        for (int r = 0; r < this.rows; r++) {
          float hold = 0;
          for (int i =0; i < this.cols; i++) {
            hold += this.data[i][r]*matrix.data[c][i];
          }
          result.data[c][r] = hold;
        }
      }
      return result;
    }
    // Note: exclude_row and exclude_col starts from 1
    Matrix subMatrix(Matrix matrix, int exclude_row, int exclude_col) {
        Matrix result = new Matrix(matrix.rows - 1, matrix.cols - 1);

        for (int row = 0, p = 0; row < matrix.rows; ++row) {
            if (row != exclude_row - 1) {
                for (int col = 0, q = 0; col < matrix.cols; ++col) {
                    if (col != exclude_col - 1) {
                        result.data[p][q] = matrix.data[row][col];

                        ++q;
                    }
                }

                ++p;
            }
        }

        return result;
    }

    float determinant() {
        if (rows != cols) {
            return -.01;
        }
        else {
            return _determinant(this);
        }
    }

    float _determinant(Matrix matrix) {
        if (matrix.cols == 1) {
            return matrix.data[0][0];
        }
        else if (matrix.cols == 2) {
            return (matrix.data[0][0] * matrix.data[1][1] -
                    matrix.data[0][1] * matrix.data[1][0]);
        }
        else {
            float result = 0.0;

            for (int col = 0; col < matrix.cols; ++col) {
                Matrix sub = subMatrix(matrix, 1, col + 1);

                result += (Math.pow(-1, 1 + col + 1) *
                           matrix.data[0][col] * _determinant(sub));
            }

            return result;
        }
    }

    Matrix inverse() {
        float det = determinant();

        if (rows != cols || det == 0.0) {
            return null;
        }
        else {
            Matrix result = new Matrix(rows, cols);

            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    Matrix sub = subMatrix(this, row + 1, col + 1);

                    result.data[col][row] = (1.0 / det *
                                             pow(-1, row + col) *
                                             _determinant(sub));
                }
            }

            return result;
        }
    }
}
