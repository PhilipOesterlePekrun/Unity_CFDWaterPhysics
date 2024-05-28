using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra.Double;
using static UnityEditor.Progress;
using static UnityEditor.Searcher.SearcherWindow.Alignment;
using System.IO;
using MathNet.Numerics;
using static System.Math;
using Unity.VisualScripting;

struct fieldExtrema
{
    public double max;
    public int iMax;
    public int jMax;
    public int nMax;
    public double min;
    public int iMin;
    public int jMin;
    public int nMin;

    public void print()
    {
        Debug.Log($"Field Extrema: max " + this.max + " @ " + this.iMax + ", " + this.jMax + ", " + this.nMax + " | min " + this.min + " @ " + this.iMin + ", " + this.jMin + ", " + this.nMin);
    }
}

// // uses chorin's projection method
public class FDM_Incomp_Inviscid_2D_chorinsProj : MonoBehaviour
{
    // external/visual
    public GameObject cameraObject;
    public GameObject backgroundObject;
    public GameObject vectorVisualsPrefab;
    public GameObject vectorVisualsParent;
    GameObject[,] vectorVisuals;
    public GameObject pVisualsPrefab; // cube or plane or something
    public GameObject pVisualsParent;
    GameObject[,] pVisuals; // pressure visuals

    double vectorScaleFactor = 4;
    double displayMultiplier = 100;
    double replaySpeed =0.1;
    string txtFileName = "Testing1";

    // physical parameters
    double[] bodyForce =new double[] { 0, -9.81 };
    const double density =998.23;
    double viscosity_mu = 1.002 * System.Math.Pow(10, -3);
    double viscosity_nu;

    // numerical parameters
    const int xCount = 40;
    const double dx = 0.05;
    double[] xBounds; // rectangle; implicitly given

    const int yCount = 40;
    const double dy = 0.05;
    double[] yBounds; //

    const int nCount = 100000; // num of time steps
    const double dt = 0.0001;
    double t_f;

    // iteration
    int iterMax = 50; // max iter steps

    // global variables
    double[,,] U; // U in i, j, n
    double[,,] V; //
    double[,,] P; //
    fieldExtrema pE;


    // // misc utilities
    void fEF(double[,,] field,ref fieldExtrema fE, int i, int j, int n) // field extrema function
    {
        if (field[i, j, n] > fE.max)
        {
            fE.max = field[i, j, n];
            fE.iMax = i;
            fE.jMax = j;
            fE.nMax = n;
        }
        if (field[i, j, n] < fE.min)
        {
            fE.min = field[i, j, n];
            fE.iMin = i;
            fE.jMin = j;
            fE.nMin = n;
        }
    }

    float function_sigmoid(float x)
    {
        return 1 / (1 + Mathf.Exp(-x));
    }
    double function_sigmoid(double x)
    {
        return 1 / (1 + Exp(-x));
    }
    double function_step(double x)
    {
        if (x < 0)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }

    double xYToAngle(double x, double y)
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
            if (y == 0)
            {
                return 0;
            }
            else
            {
                return 180 / PI * Atan(y / x);
            }
        }
        else // == else if(x < 0)
        {
            if (y == 0)
            {
                return 180;
            }
            else
            {
                return 180 + 180 / PI*Atan(y / x);
            }
        }
    }

    // // scenario functions (e.g. BC velocities in time and along boundary or body force in time and position)
    double U_onUpperBoundary(double x,double t)
    {
        double frequ =0.5; // Hz
        return 0*Max(0.7*function_step(t-0.8),-2*Pow(t-0.8,2)+0.7);//0.4f *(2*function_sigmoid(t)-1);//System.Math.Sin(2 * System.Math.PI * frequ*t);
    }
    double U_onLeftBoundary(double y,double t)
    {
        return 0*Max(0.7 * function_step(t - 0.8), -2 * Pow(t - 0.8, 2) + 0.7);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Start()
    {
        // grid
        xBounds = new double[] { 0, (xCount - 1) * dx };
        yBounds = new double[] { 0, (yCount - 1) * dy };
        t_f = (nCount - 1) * dt;
        Debug.Log("x length,y length, t_f: " + (xBounds[1] - xBounds[0]) + "m , " + (yBounds[1] - yBounds[0]) + "m , " + t_f + "s");

        // material properties and constants
        viscosity_nu = density * viscosity_mu;

        // vectors and visuals
        U = new double[xCount, yCount, nCount]; // x-component of velocity
        V = new double[xCount, yCount, nCount]; // y-component of velocity
        P = new double[xCount, yCount, nCount]; // pressure

        pVisuals = new GameObject[xCount, yCount];
        vectorVisuals = new GameObject[xCount, yCount];

        // // ICs: uniform =0
        double p_0 = 0; // reference pressure @ y=0m
        for (int i = 0; i < xCount; i++)
        {
            for (int j = 0; j < yCount; j++)
            {
                U[i, j, 0] = 0;//Max(0, -Pow(i -10, 2)-Pow(j - yCount / 2, 2) + 0.8);
                V[i, j, 0] = 0;
                P[i, j, 0] = 0;
                //P[i, j, 0] =p_0 + density * bodyForce[1] * (j * dy);

                // visuals:
                vectorVisuals[i, j] = Object.Instantiate(vectorVisualsPrefab, new Vector3((float)displayMultiplier * (float)(i * dx), (float)displayMultiplier * (float)(j * dy), 0), Quaternion.identity, vectorVisualsParent.transform);
                vectorVisuals[i, j].transform.position += new Vector3(0, 0, -1);

                pVisuals[i, j] = Object.Instantiate(pVisualsPrefab, new Vector3((float)displayMultiplier * (float)(i * dx), (float)displayMultiplier * (float)(j * dy), 0), Quaternion.identity, pVisualsParent.transform);
                pVisuals[i, j].transform.localScale = new Vector3((float)displayMultiplier*(float)dx, (float)displayMultiplier*(float)dy,1);
            }
        }
        // pressure ICBCs
        for (int i = 1; i < xCount - 1; i++) // top and bottom impermeable
        {
            P[i, 0,0] = P[i, 1,0];
            P[i, yCount - 1,0] = P[i, yCount - 2,0];
        }
        for (int j = 1; j < yCount - 1; j++) // right and left impermeable
        {
            P[0, j,0] = P[1, j,0];
            P[xCount - 1, j,0] = P[xCount - 2, j,0];
        }

        // // time loop
        double dx2 = Pow(dx, 2); // calculated only once for efficiency
        double dy2 = Pow(dy, 2); //
        for (int n = 1; n < nCount; n++)
        {
            // // 1. intermediate velocity
            double[,] U_star = new double[xCount, yCount];
            double[,] V_star = new double[xCount, yCount];

            //BCs
            for (int i = 1; i < xCount - 1; i++) // top and bottom
            {
                U_star[i, yCount - 1] = U_onUpperBoundary(xBounds[0] + i * dx, n * dt);
                V_star[i, yCount - 1] = 0;
                U_star[i, 0] = 0;
                V_star[i, 0] = 0;
            }
            for (int j = 1; j < yCount - 1; j++) // right and left
            {
                U_star[xCount - 1, j] = U_onLeftBoundary(xBounds[0] + j* dy, n * dt);
                V_star[xCount - 1, j] = 0;
                U_star[0, j] = U_onLeftBoundary(xBounds[0] + j * dy, n * dt);
                V_star[0, j] = 0;
            }
            for (int i =1; i < xCount-1; i++)
            {
                for (int j =1; j < yCount-1; j++)
                {
                    U_star[i, j] = U[i, j, n - 1] + dt * (-(U[i, j, n - 1] * (U[i + 1, j, n - 1] - U[i - 1, j, n - 1]) / (2 * dx) + V[i, j, n - 1] * (U[i, j + 1, n - 1] - U[i, j - 1, n - 1]) / (2 * dy)) + bodyForce[0]);
                    V_star[i, j] = V[i, j, n - 1] + dt * (-(U[i, j, n - 1] * (V[i + 1, j, n - 1] - V[i - 1, j, n - 1]) / (2 * dx) + V[i, j, n - 1] * (V[i, j + 1, n - 1] - V[i, j - 1, n - 1]) / (2 * dy)) + bodyForce[1]);
                }
            }

            // // 2. iterate poisson eq
            // initial pressure guess = pressure from time step n-1
            for (int i = 0; i < xCount; i++)
            {
                for (int j = 0; j < yCount; j++)
                {
                    P[i, j, n] = P[i, j, n - 1];
                }
            }
            // iterate
            double[,] pNew = new double[xCount, yCount];
            for (int iter = 0; iter < iterMax; iter++)
            {
                for (int i = 1; i < xCount - 1; i++)
                {
                    for (int j = 1; j < yCount - 1; j++)
                    {
                        pNew[i, j] = 1 / (2* (dx2 + dy2)) * ((P[i + 1, j, n] + P[i - 1, j, n]) * dy2 + (P[i, j + 1, n] + P[i, j - 1, n]) * dx2 -
                            dx2 * dy2 * density / dt * ((U_star[i + 1, j] - U_star[i - 1, j]) / (2 * dx) + (V_star[i, j + 1] - V_star[i, j - 1]) / (2 * dy)));
                    }
                }

                // pressure BCs
                for (int i = 1; i < xCount - 1; i++) // top and bottom impermeable
                {
                    pNew[i, 0] = pNew[i, 1];
                    pNew[i, yCount - 1] = pNew[i, yCount - 2];
                }
                for (int j = 1; j < yCount - 1; j++) // right and left impermeable
                {
                    pNew[0, j] = pNew[1, j];
                    pNew[xCount - 1, j] = pNew[xCount - 2, j];
                }

                for (int i = 0; i < xCount; i++)
                {
                    for (int j = 0; j < yCount; j++)
                    {
                        P[i, j, n] = pNew[i, j];
                    }
                }
            }

            // // 3. get u^{n+1}
            //BCs
            for (int i = 1; i < xCount - 1; i++) // top and bottom
            {
                U[i, yCount - 1, n] = U_onUpperBoundary(xBounds[0] + i * dx, n * dt);
                V[i, yCount - 1, n] = 0;
                U[i, 0, n] = 0;
                V[i, 0, n] = 0;
            }
            for (int j = 1; j < yCount - 1; j++) // right and left
            {
                U[xCount - 1, j, n] = U_onLeftBoundary(xBounds[0] + j * dy, n * dt);
                V[xCount - 1, j, n] = 0;
                U[0, j, n] = U_onLeftBoundary(xBounds[0] + j * dy, n * dt);
                V[0, j, n] = 0;
            }
            for (int i = 1; i < xCount-1; i++)
            {
                for (int j = 1; j < yCount-1; j++)
                {
                    U[i, j, n] = U_star[i, j] - dt / density * (P[i + 1, j, n] - P[i - 1, j, n]) / (2 * dx);
                    V[i, j, n] = V_star[i, j] - dt / density * (P[i, j + 1, n] - P[i, j - 1, n]) / (2 * dy);
                }
            }
        }

        // get min and max values
        fieldExtrema uE = new fieldExtrema();
        fieldExtrema vE = new fieldExtrema();
        pE = new fieldExtrema();
        uE.max = double.MinValue;
        uE.min = double.MaxValue;
        vE.max = double.MinValue;
        vE.min = double.MaxValue;
        pE.max = double.MinValue;
        pE.min = double.MaxValue;
        for (int n = 0; n < nCount; n++)
        {
            for (int i = 0; i < xCount; i++)
            {
                for (int j = 0; j < yCount; j++)
                {
                    fEF(U, ref uE, i, j, n);
                    fEF(V, ref vE, i, j, n);
                    if (!((i == 0 || i == xCount - 1) && (j == 0 && j == yCount - 1))) // leave the corners out of it
                    {
                        fEF(P, ref pE, i, j, n);
                    }
                }
            }
        }
        Debug.Log("uE "); uE.print();
        Debug.Log("vE "); vE.print();
        Debug.Log("pE "); pE.print();
        Debug.Log("UONUPPER " +U_onUpperBoundary(10*dx,1.2));
    }


    // Update is called once per frame
    void Update()
    {
        if (Time.time > 2)
        {
            int n = 0;
            n = ((int)((Time.time - 2) / dt * replaySpeed)) %(nCount-1);
            Debug.Log("t= "+n * dt);
            Debug.Log("U @ i=10, j=yCount-1: "+(float)U[10,yCount-1, n]);
            for (int i = 0; i < xCount; i++)
            {
                for (int j = 0; j < yCount; j++)
                {
                    GameObject pV = pVisuals[i, j];
                    MeshRenderer[] mR = pV.GetComponentsInChildren<MeshRenderer>();
                    for (int k = 0; k < mR.Length; k++)
                    {
                        mR[k].material.color = new Color((float)(0.25+0.5*(P[i, j, n] - pE.min) / (pE.max - pE.min)),0,0);
                    }

                    GameObject vV = vectorVisuals[i,j];

                    double magnitude = Max(-20,Min(20,Sqrt(Pow(U[i, j, n], 2) +Pow(V[i, j, n], 2))));

                    vV.transform.localScale =new Vector3((float)vectorScaleFactor*(float)magnitude, 1, 1);
                    vV.transform.rotation = Quaternion.Euler(0, 0,(float)xYToAngle(U[i, j, n], V[i,j,n]));

                    /*
                    MeshRenderer[] mR = vV.GetComponentsInChildren<MeshRenderer>();
                    for (int i = 0; i < mR.Length; i++)
                    {
                        mR[i].material.color = new Color(function_sigmoid(magnitude - 2), function_sigmoid(magnitude - 2), function_sigmoid(magnitude - 2));
                    }*/
                }
            }
        }
    }
}
