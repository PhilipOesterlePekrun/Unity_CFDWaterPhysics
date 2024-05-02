using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra.Double;
using static UnityEditor.Progress;
using static UnityEditor.Searcher.SearcherWindow.Alignment;
using System.IO;
using MathNet.Numerics;
using static System.Math;

public class FDM_Incomp_Viscous_StaggeredVelocity : MonoBehaviour
{
    public GameObject cameraObject;
    public GameObject backgroundObject;
    public GameObject vectorVisualPrefab;
    public GameObject vectorVisualsParent;
    GameObject[] vectorVisuals;
    int displayMultiplier = 150;
    string txtFileName = "Testing1";
    string iterExampleLogFileName = "iterLogExample";
    int ILE_I = 10;
    int ILE_J = 10;
    int ILE_N = 10;
    float vectorScaleFactor = 1;
    float replaySpeed = 0.1f;

    DenseVector bodyForce = new DenseVector(new double[] { 0, -9.81 });
    double density = 998.23;
    double viscosity_mu = 1.002 * System.Math.Pow(10, -3);
    double viscosity_nu;

    // numerical parameters
    int xCount = 20; // num of nodes; implicitly given
    double dx = 0.05;
    double[] xBounds; // rectangle

    int yCount = 20; //
    double dy = 0.05;
    double[] yBounds; //

    int nCount = 20; // num of time steps
    double dt = 0.01;
    double t_f;

    int iter_maxSteps = 15; // maximal number of iteration steps for implicit equation solving
    double iter_tuningCoefficient = 0.1; // [Pa/step]
    double iter_minError = 0.05; // if lower than this error, no further steps will be taken [Pa]

    DenseVector[] AInTimeGlobal;

    // // misc utilities
    float function_sigmoid(float x)
    {
        return 1 / (1 + Mathf.Exp(-x));
    }
    double function_sigmoid(double x)
    {
        return (double)function_sigmoid((float)x);
    }

    // // transformations
    int vectToLin(int i, int j) // lexicographic and starts from k=1
    {
        return ((j - 1) * xCount) + i;
    }
    int[] linToVect(int k) // lexicographic as above, starts from i=1 and j=1
    {
        return new int[] { (k - 1) % xCount + 1, 1 + Mathf.FloorToInt((k - 1) / xCount) };
    }

    float xYToAngle(double x, double y)
    {
        if (x == 0)
        {
            if (y > 0)
            {
                return 90;
            }
            if (y < 0)
            {
                return 270;
            }
            else
            {
                return 0;
            }
        }
        if (x > 0) // == else if
        {
            return 180 / Mathf.PI * Mathf.Atan((float)y / ((float)x));
        }
        else // == else if(x < 0)
        {
            return 180 / Mathf.PI * (Mathf.PI - Mathf.Atan((float)y / ((float)x)));
        }
    }

    void Start()
    {
        // grid
        xBounds = new double[] { 0, xCount * dx };
        yBounds = new double[] { 0, yCount * dy };
        t_f = nCount * dt;
        Debug.Log("x length,y length, t_f: " + (xBounds[1] - xBounds[0]) + "m , " + (yBounds[1] - yBounds[0]) + "m , " + t_f + "s");

        // material properties and constants
        viscosity_nu = density * viscosity_mu;

        // vectors and visuals
        DenseVector A_0 = new DenseVector(3 * vectToLin(xCount, yCount));
        vectorVisuals = new GameObject[3 * vectToLin(xCount, yCount)];
        // ICs: hydrostatic equilibrium
        double p_0 = 0; // reference pressure @ y=0m
        for (int k = 1; k <= vectToLin(xCount, yCount); k++)
        {
            int[] lV = linToVect(k);
            if (lV[0] == 1 || lV[0] == xCount || lV[1] == 1 || lV[1] == yCount) // BC
            {
                A_0.Values[3 * (k - 1)] = 0;
                A_0.Values[3 * (k - 1) + 1] = 0;
            }
            //else if // set some initial velocities in some area (smooth, continuous)
            else
            {
                A_0.Values[3 * (k - 1)] = (double)(Mathf.Max(0, 0.2f - 0.02f * (0.2f * Mathf.Pow((lV[0] - xCount / 2), 2) + Mathf.Pow((lV[1] - yCount / 2), 2))));
                A_0.Values[3 * (k - 1) + 1] = (double)(Mathf.Max(0, 0.2f - 0.05f * (Mathf.Pow((lV[0] - xCount / 2), 2) + Mathf.Pow((lV[1] - yCount / 2), 2))));
            }
            A_0.Values[(k - 1) + 2] = p_0 + density * bodyForce.Values[1] * linToVect(k)[1];

            // vector visuals:
            vectorVisuals[k] = Object.Instantiate(vectorVisualPrefab, new Vector3(displayMultiplier * (float)((lV[0] - 1) * dx), displayMultiplier * (float)((lV[1] - 1) * dy), 0), Quaternion.identity, vectorVisualsParent.transform);
        }

        DenseVector[] AInTime = new DenseVector[nCount];
        AInTime[0] = A_0;
        double[] ILE = new double[15];
        for (int n = 1; n < nCount; n++) // t = n * dt // start from second time step which is n=1
        {
            DenseVector ALast = AInTime[n - 1];
            // iteration with maximum step count iter_maxSteps // if the divergence is positive, the pressure should be increased to compensate
            int iter_step = 0;
            DenseVector ATemp = new DenseVector(3 * vectToLin(xCount, yCount)); // this time step
            for (int k = 1; k <= vectToLin(xCount, yCount); k++)
            {
                int[] lV = linToVect(k);
                if (lV[0] == 1 || lV[0] == xCount || lV[1] == 1 || lV[1] == yCount) // BC
                {
                    ATemp.Values[3 * (k - 1)] = 0;
                    ATemp.Values[3 * (k - 1) + 1] = 0;
                }
                else
                {
                    double iter_pressure = ALast.Values[3 * (k - 1) + 2];
                    double divUV = iter_minError + 1; // initial value will not be used but should not stop loop, i.e. can't be ~0

                    bool bool1 = false;
                    if (ILE_I == lV[0] && ILE_J == lV[1] && ILE_N == n)
                    {
                        bool1 = true;
                    }
                    while (iter_step < iter_maxSteps && divUV > iter_minError)
                    {
                        double eqU_1 = -(ALast.Values[3 * (k - 1)] * (ALast.Values[3 * (vectToLin(lV[0] + 1, lV[1]) - 1)] - ALast.Values[3 * (vectToLin(lV[0] - 1, lV[1]) - 1)]) / (2 * dx) +
                        ALast.Values[3 * (k - 1) + 1] * (ALast.Values[3 * (vectToLin(lV[0], lV[1] + 1) - 1)] - ALast.Values[3 * (vectToLin(lV[0], lV[1] - 1) - 1)]) / (2 * dy));
                        double eqV_1 = -(ALast.Values[3 * (k - 1) + 1] * (ALast.Values[3 * (vectToLin(lV[0] + 1, lV[1]) - 1) + 1] - ALast.Values[3 * (vectToLin(lV[0] - 1, lV[1]) - 1) + 1]) / (2 * dx) +
                            ALast.Values[3 * (k - 1) + 1] * (ALast.Values[3 * (vectToLin(lV[0], lV[1] + 1) - 1) + 1] - ALast.Values[3 * (vectToLin(lV[0], lV[1] - 1) - 1) + 1]) / (2 * dy));
                        double eqU_2 = viscosity_nu * ((ALast.Values[3 * (vectToLin(lV[0] + 1, lV[1]) - 1)] - 2 * ALast.Values[3 * (k - 1)] + ALast.Values[3 * (vectToLin(lV[0] - 1, lV[1]) - 1)]) / System.Math.Pow(dx, 2) +
                            (ALast.Values[3 * (vectToLin(lV[0], lV[1] + 1) - 1)] - 2 * ALast.Values[3 * (k - 1)] + ALast.Values[3 * (vectToLin(lV[0], lV[1] - 1) - 1)]) / System.Math.Pow(dy, 2));
                        double eqV_2 = viscosity_nu * ((ALast.Values[3 * (vectToLin(lV[0] + 1, lV[1]) - 1) + 1] - 2 * ALast.Values[3 * (k - 1) + 1] + ALast.Values[3 * (vectToLin(lV[0] - 1, lV[1]) - 1) + 1]) / System.Math.Pow(dx, 2) +
                            (ALast.Values[3 * (vectToLin(lV[0], lV[1] + 1) - 1) + 1] - 2 * ALast.Values[3 * (k - 1) + 1] + ALast.Values[3 * (vectToLin(lV[0], lV[1] - 1) - 1)] + 1) / System.Math.Pow(dy, 2));
                        double eqU_3 = -1 / density * (ALast.Values[3 * (vectToLin(lV[0] + 1, lV[1]) - 1) + 2] - ALast.Values[3 * (vectToLin(lV[0] - 1, lV[1]) - 1) + 2]) / (2 * dx);
                        double eqV_3 = -1 / density * (ALast.Values[3 * (vectToLin(lV[0], lV[1] + 1) - 1) + 2] - ALast.Values[3 * (vectToLin(lV[0], lV[1] - 1) - 1) + 2]) / (2 * dx);

                        ATemp.Values[3 * (k - 1)] = ALast.Values[3 * (k - 1)] + dt * eqU_1 + eqU_2 + eqU_3;
                        ATemp.Values[3 * (k - 1) + 1] = ALast.Values[3 * (k - 1)] + dt * eqV_1 + eqV_2 + eqV_3;

                        divUV = (ALast.Values[3 * (vectToLin(lV[0] + 1, lV[1]) - 1)] - ALast.Values[3 * (vectToLin(lV[0] - 1, lV[1]) - 1)]) / (2 * dx) +
                        (ALast.Values[3 * (vectToLin(lV[0], lV[1] + 1) - 1)] - ALast.Values[3 * (vectToLin(lV[0], lV[1] - 1) - 1)]) / (2 * dy);
                        if (bool1)
                        {
                            ILE[iter_step] = iter_pressure;
                        }
                        iter_pressure -= iter_tuningCoefficient * divUV;

                        iter_step++;
                    }
                }
            }
            AInTime[n] = ATemp;
        }
        AInTimeGlobal = AInTime;

        // // write to file
        string filePath = Application.dataPath + @"\" + txtFileName + ".txt";
        Debug.Log(filePath);
        using (StreamWriter writer = new StreamWriter(filePath))
        {
            for (int n = 0; n < AInTimeGlobal.Length; n++)
            {
                writer.WriteLine("");
                writer.WriteLine("TIME STEP N = " + n);
                for (int k = 1; k <= vectToLin(xCount, yCount); k++)
                {
                    int[] lV = linToVect(k);
                    writer.WriteLine("k = " + k + " | i, j = " + lV[0] + ", " + lV[1] + " | U, V = " + AInTimeGlobal[n].Values[2 * (k - 1)] + ", " + AInTimeGlobal[n].Values[2 * (k - 1) + 1]);
                }
            }
        }
        // // write iter log example
        string filePathILE = Application.dataPath + @"\" + iterExampleLogFileName + ".txt";
        Debug.Log(filePath);
        using (StreamWriter writer = new StreamWriter(filePathILE))
        {
            writer.WriteLine("Iter Log Example for " + "i= " + ILE_I + ", j= " + ILE_J + ", n= " + ILE_N);
            for (int i = 0; i < ILE.Length; i++)
            {
                writer.WriteLine(ILE[i]);
            }
        }

        for (int k = 1; k <= vectToLin(xCount, yCount); k++)
        {
            GameObject vV = vectorVisuals[k];
        }
    }
    void Update()
    {
        if (Time.time > 2)
        {
            int n = 0;
            n = ((int)((Time.time - 2) / dt * replaySpeed)) % nCount;
            Debug.Log(n * dt);

            for (int k = 1; k <= vectToLin(xCount, yCount); k++)
            {
                GameObject vV = vectorVisuals[k];

                float magnitude = Mathf.Max(-10, Mathf.Min(10, Mathf.Sqrt(Mathf.Pow((float)AInTimeGlobal[n].Values[2 * (k - 1)], 2) + Mathf.Pow((float)AInTimeGlobal[n].Values[2 * (k - 1) + 1], 2))));

                vV.transform.localScale = vectorScaleFactor * new Vector3(magnitude, 1, 1);
                vV.transform.rotation = Quaternion.Euler(0, 0, -xYToAngle(AInTimeGlobal[n].Values[3 * (k - 1)], AInTimeGlobal[n].Values[3 * (k - 1) + 1]));

                MeshRenderer[] mR = vV.GetComponentsInChildren<MeshRenderer>();
                for (int i = 0; i < mR.Length; i++)
                {
                    mR[i].material.color = new Color(function_sigmoid(magnitude - 2), function_sigmoid(magnitude - 2), function_sigmoid(magnitude - 2));
                }
            }
        }
    }
}