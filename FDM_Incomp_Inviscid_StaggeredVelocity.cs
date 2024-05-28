using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra.Double;
using static UnityEditor.Progress;
using static UnityEditor.Searcher.SearcherWindow.Alignment;
using System.IO;
using MathNet.Numerics;
using static System.Math;

public class FDM_Incomp_Inviscid_2D_StaggeredVelocity: MonoBehaviour
{
    public GameObject cameraObject;
    public GameObject backgroundObject;
    public GameObject vectorVisualPrefab;
    public GameObject vectorVisualsParent;
    public GameObject cubePrimitivePrefab;
    public GameObject pVisualsParent;
    //GameObject[] vectorVisualsU; // velocity; staggered; arrows
    //GameObject[] vectorVisualsV; // velocity; staggered; arrows
    GameObject[] vectorVisuals; // actually I will just have the averaged total velocity arrow at the pressure nodes as the velocity visual
    GameObject[] pVisuals; // pressure; standard nodes; colored rectangles
    int displayMultiplier = 150;
    string txtFileName = "Testing1";
    string iterExampleLogFileName = "iterLogExample";
    int ILE_I = 10; // Iteratio
    int ILE_J = 10;
    int ILE_N = 10;
    float vectorScaleFactor = 1;
    float replaySpeed = 0.1f;

    DenseVector bodyForce = new DenseVector(new double[] { 0, -9.81 });
    const double density = 998.23;
    double viscosity_mu = 1.002 * System.Math.Pow(10, -3);
    double viscosity_nu;

    // numerical parameters
    // NOTE: I call the pressure nodes the standard nodes, or just nodes, while staggered nodes are specifically denoted as such, usually by variables and methods ending in S
    const int xCount = 20; // num of pressure nodes
    const double dx = 0.05;
    double[] xBounds; // rectangle; implicitly given

    const int yCount = 20; //
    const double dy = 0.05;
    double[] yBounds; //

    const int nCount = 20; // num of time steps
    const double dt = 0.01;
    double t_f;

    const int iter_maxSteps = 15; // maximal number of iteration steps for implicit equation solving
    const double iter_tuningCoefficient = 0.1; // [Pa/step]
    const double iter_minError = 0.05; // if lower than this error, no further steps will be taken [Pa]

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
        return ((j - 1) * xCount) +(i-1);
    }
    int vectToLinSU(int i, int j) // staggered U variant
    {
        return ((j - 1) * xCount) +(i);
    }
    int vectToLinSV(int i, int j) // staggered V variant
    {
        return ((j) * xCount) +(i-1);
    }
    int[] linToVect(int k) // lexicographic as above, starts from i=1 and j=1
    {
        return new int[] {1+(k - 1) % xCount, 1 + Mathf.FloorToInt((k - 1) / xCount) };
    }
    int[] linToVectSU(int k) // staggered U variant
    {
        return new int[] {1+(k - 1) %(xCount), 1 + Mathf.FloorToInt((k - 1) / xCount) };
    }
    int[] linToVectSV(int k) // staggered V variant
    {
        return new int[] {1+(k - 1) %(xCount), 1 + Mathf.FloorToInt((k - 1) / xCount) };
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

    // // scenario functions (e.g. BC velocities in time and along boundary or body force in time and position)
    double U_onUpperBoundary(float x,float t)
    {
        double frequ = 0.2; // Hz
        return 0.1f *System.Math.Sin(2 * System.Math.PI * frequ);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void Start()
    {
        //var mat=MathNet.Numerics.LinearAlgebra.CreateMatrix.Dense<float>(5,5);
        //var inv=mat.Inverse();
        // grid
        xBounds = new double[] { 0,xCount* dx };
        yBounds = new double[] { 0,yCount* dy };
        t_f =(nCount+1)* dt;
        Debug.Log("x length,y length, t_f: " + (xBounds[1] - xBounds[0]) + "m , " + (yBounds[1] - yBounds[0]) + "m , " + t_f + "s");

        // material properties and constants
        viscosity_nu = density * viscosity_mu;

        // vectors and visuals
        DenseVector U_0= new DenseVector(vectToLinSU(xCount, yCount-1)); // x-component of velocity at staggered nodes
        DenseVector V_0= new DenseVector(vectToLinSV(xCount-1, yCount)); // y-component " "
        DenseVector P_0= new DenseVector(vectToLin(xCount-1, yCount-1)); // pressure at pressure/standard nodes
        pVisuals= new GameObject[vectToLin(xCount-1, yCount-1)];
        vectorVisuals = new GameObject[vectToLin(xCount - 1, yCount - 1)];
        // ICs: hydrostatic equilibrium
        double p_0 = 0; // reference pressure @ y=0m
        for (int k = 1; k <= vectToLin(xCount, yCount); k++)
        {
            int[] lV = linToVect(k);
            
            P_0.Values[k]= p_0 + density * bodyForce.Values[1] * linToVect(k)[1];

            // visuals:
            pVisuals[k] = Object.Instantiate(cubePrimitivePrefab, new Vector3(displayMultiplier * (float)((lV[0] - 1) * dx), displayMultiplier * (float)((lV[1] - 1) * dy), 0), Quaternion.identity, vectorVisualsParent.transform);
            vectorVisuals[k] = Object.Instantiate(vectorVisualPrefab, new Vector3(displayMultiplier * (float)((lV[0] - 1) * dx), displayMultiplier * (float)((lV[1] - 1) * dy), 0), Quaternion.identity, vectorVisualsParent.transform);
        } ////LEFT OFF HERE
        for(int k = 1; k <= vectToLinSU(xCount, yCount+1); k++) // SU
        {
            int[] lV = linToVectSU(k);
            if (lV[0] == 1 || lV[0] == xCount || lV[1] == 1) // on boundary except upper
            {
                U_0.Values[k] = 0;
            }
            else if(lV[1] == yCount + 1) // on upper boundary
            {
                U_0.Values[k] = U_onUpperBoundary(lV[0], 0);
            }
            //else if // set some initial velocities in some area (smooth, continuous)
            else
            {
                U_0.Values[k] = 0;
                //U_0[k]= (double)(Mathf.Max(0, 0.2f - 0.02f * (0.2f * Mathf.Pow((lV[0] - xCount / 2), 2) + Mathf.Pow((lV[1] - yCount / 2), 2)))); // ACTUALLY I will just set some velocity boundary conditions instead, for now
            }
            //vectorVisualsU[k] = Object.Instantiate(vectorVisualPrefab, new Vector3(displayMultiplier * (float)((lV[0] - 1+0.5f) * dx), displayMultiplier * (float)((lV[1] - 1) * dy), 0), Quaternion.identity, vectorVisualsParent.transform);
        } ////LEFT OFF HERE
        for (int k = 1; k <= vectToLinSV(xCount+1, yCount); k++) // SV
        {
            int[] lV = linToVectSV(k);
            if (lV[0] == 1 || lV[0] == xCount+1 || lV[1] == 1 || lV[1] == yCount) // on boundary
            {
                V_0.Values[k] = 0;
            }
            //else if // set some initial velocities in some area (smooth, continuous)
            else
            {
                V_0.Values[k] = 0;
                //V_0[k] = (double)(Mathf.Max(0, 0.2f - 0.05f * (Mathf.Pow((lV[0] - xCount / 2), 2) + Mathf.Pow((lV[1] - yCount / 2), 2)))); //
            }
            //vectorVisualsV[k] = Object.Instantiate(vectorVisualPrefab, new Vector3(displayMultiplier * (float)((lV[0] - 1) * dx), displayMultiplier * (float)((lV[1] - 1+0.5f) * dy), 0), Quaternion.identity, vectorVisualsParent.transform);
        } ////LEFT OFF HERE

        // assemble A_0
        DenseVector A_0 = new DenseVector(3*Mathf.Max(U_0.Count, V_0.Count, P_0.Count)); // will not be entirely filled because there are fewer pressure than velocity nodes, and SU and SV node count might also differ
        for(int k =0; k < A_0.Count; k+=3)
        {
            A_0.Values[k] = U_0.Values[k / 3];
            A_0.Values[k+1] = U_0.Values[k / 3];
            A_0.Values[k+2] = U_0.Values[k / 3];
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

/*
TODO:
- newton-raphson iteration
- convert all float operations and variables to double (use System.Math, which is included thus the call of e.g. double Sin(double x) is a call to System.Math.Sin()
*/