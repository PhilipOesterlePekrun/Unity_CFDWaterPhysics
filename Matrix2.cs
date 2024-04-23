public class Matrix2 // 2nd order matrix, though vectors are essentially a subset of this
{
    float[,] m;
    int rowCount = 1;
    int colCount = 1;
    public Matrix2(int rowCount, int colCount) // true 2nd order matrix, eg a 2x3 matrix is Matrix2(2,3)
    {
        m= new float[rowCount, colCount];
        this.rowCount = rowCount;
        this.colCount = colCount;
    }
    public Matrix2(int rowCount) // vector
    {
        m= new float[rowCount,1];
        this.rowCount = rowCount;
        this.colCount = 1;
    }
    float getEle(int pos1, int pos2)
    {
        return m[pos1, pos2];
    }
    float getEle(int pos)
    {
        return m[pos,1];
    }
    // Start is called before the first frame update
    void changeEle(int pos1, int pos2, float ele)
    {
        m[pos1, pos2] = ele;
    }
    void changeEle(int pos, float ele)
    {
        m[pos, 1] = ele;
    }
    void transpose()
    {
        float[,] newM= new float[colCount, rowCount];
        for(int i = 0; i < rowCount;i++)
        {
            for(int j = 0; j < colCount; j++)
            {
                newM[j, i] = m[i, j];
            }
        }
        m = newM;
    }
    float determinant()
    {
        if (rowCount != colCount)
        {
            //Debug.Log("__ERROR: Invalid Determinant");
            return 0;
        }
        else if (rowCount == 1)
        {
            return m[0, 0];
        }
        else if (rowCount == 2)
        {
            return m[0, 0] * m[1, 1] - m[1, 2] * m[2, 1];
        }
        else if (rowCount == 3)
        {
            return m[0, 0] * m[1, 1] * m[2, 2] + m[0,1 ] * m[1, 2] * m[2, 0] + m[0, 2] * m[1, 0] * m[2,1] - (m[2, 0] * m[1, 1] * m[0, 2] + m[2, 1] * m[1, 2] * m[0, 0] + m[2, 2] * m[1, 0] * m[0,1]);
        }
        else
        {
            return 0;/////////////////////////////left off
        }
    }
}
