using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra.Double;
using static UnityEditor.Progress;
using static UnityEditor.Searcher.SearcherWindow.Alignment;
using System.IO;
using MathNet.Numerics;

//// DISCONTINUED BEFORE DOING PRESSURE
public class FDM_Incomp_Inviscid_Unsteady_2D_noPressure : MonoBehaviour
{
    public GameObject cameraObject;
    public GameObject backgroundObject;
    public GameObject vectorVisualPrefab;
    public GameObject vectorVisualsParent;
    GameObject[] vectorVisuals;
    int displayMultiplier = 150;
    string txtFileName = "Testing1";
    public Material white;
    float vectorScaleFactor=1;
    float replaySpeed=0.5f;

    DenseVector bodyForce = new DenseVector(new double[]{0, -9.81});
    double density = 998.23;
    double dx = 0.05;
    double[] xBounds = { 0,4}; // rectangle
    int xCount; // num of nodes; implicitly given
    double dy = 0.05;
    double[] yBounds = { 0,2}; //
    int yCount; //
    double dt = 0.005;
    double t_f =3;
    int nCount; // num of time steps

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
        return ((j-1)*xCount)+ i;
    }
    int[] linToVect(int k) // lexicographic as above, starts from i=1 and j=1
    {
        return new int[] {(k-1)%xCount+1,1+Mathf.FloorToInt((k-1)/xCount)};
    }

    float xYToAngle(double x,double y)
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
            return 180/Mathf.PI*Mathf.Atan((float)y /((float)x));
        }
        else // == else if(x < 0)
        {
            return 180/Mathf.PI * (Mathf.PI - Mathf.Atan((float)y / ((float)x)));
        }
    }

    void Start()
    {
        xCount = (int)((xBounds[1] - xBounds[0]) / dx)+1;
        yCount = (int)((yBounds[1] - yBounds[0]) / dy)+1;
        nCount = (int)(t_f / dt) + 1;
        Debug.Log("xCount, yCount, nCount: " + xCount + ", " + yCount + ", " + nCount);

        DenseVector A_0 = new DenseVector(3* vectToLin(xCount, yCount));
        vectorVisuals = new GameObject[3* vectToLin(xCount, yCount)];
        // ICs: hydrostatic equilibrium
        double p_0 = 0; // reference pressure @ y=0m
        for(int k = 1; k <= vectToLin(xCount, yCount); k++)
        {
            int[] lV = linToVect(k);
            if (lV[0] == 1 || lV[0] == xCount || lV[1] == 1 || lV[1] == yCount) // BC
            {
                A_0.Values[3* (k - 1)] = 0;
                A_0.Values[3* (k - 1) + 1] = 0;
            }
            //else if // set some initial velocities in some area (smooth, continuous)
            else
            {
                A_0.Values[3* (k - 1)] =(double)(Mathf.Max(0,0.2f-0.02f*(0.2f*Mathf.Pow((lV[0]-xCount/2),2)+Mathf.Pow((lV[1]-yCount/2),2))));
                A_0.Values[3* (k - 1) + 1] =(double)(Mathf.Max(0,0.2f-0.05f * (Mathf.Pow((lV[0] - xCount / 2), 2) + Mathf.Pow((lV[1] - yCount / 2), 2))));
            }
            A_0.Values[(k-1)+2] = p_0 + density*bodyForce.Values[1]*linToVect(k)[1];

            // vector visuals:
            vectorVisuals[k] = Object.Instantiate(vectorVisualPrefab, new Vector3(displayMultiplier*(float)((lV[0]-1)*dx), displayMultiplier*(float)((lV[1]-1)*dy), 0),Quaternion.identity,vectorVisualsParent.transform);
        }

        DenseVector[] AInTime = new DenseVector[nCount];
        AInTime[0] = A_0;
        for (int n = 1; n <nCount; n++) // t = n * dt // start from second time step which is n=1
        {
            DenseVector ATemp = new DenseVector(3* vectToLin(xCount, yCount)); // this time step
            for(int k = 1; k <= vectToLin(xCount, yCount); k++)
            {
                int[] lV = linToVect(k);
                if (lV[0] == 1 || lV[0] == xCount || lV[1] == 1 || lV[1] == yCount) // BC
                {
                    ATemp.Values[3* (k - 1)] = 0;
                    ATemp.Values[3* (k - 1) + 1] = 0;
                }
                else
                {
                    //with pressure term
                    DenseVector ALast = AInTime[n - 1];
                    ATemp.Values[3* (k - 1)] = ALast.Values[3* (k - 1)] + dt *
                        (-(ALast.Values[3* (k - 1)] * (ALast.Values[3* (vectToLin(lV[0] + 1, lV[1]) - 1)] - ALast.Values[3* (vectToLin(lV[0] - 1, lV[1]) - 1)]) / (2 * dx) +
                        ALast.Values[3* (k - 1) + 1] * (ALast.Values[3* (vectToLin(lV[0], lV[1] + 1) - 1)] - ALast.Values[3* (vectToLin(lV[0], lV[1] - 1) - 1)]) / (2 * dy)));
                    ATemp.Values[3* (k - 1) + 1] = ALast.Values[3* (k - 1) + 1] + dt *
                        (-(ALast.Values[3* (k - 1)] * (ALast.Values[3* (vectToLin(lV[0] + 1, lV[1]) - 1) + 1] - ALast.Values[3* (vectToLin(lV[0] - 1, lV[1]) - 1) + 1]) / (2 * dx) +
                        ALast.Values[3* (k - 1) + 1] * (ALast.Values[3* (vectToLin(lV[0], lV[1] + 1) - 1) + 1] - ALast.Values[3* (vectToLin(lV[0], lV[1] - 1) - 1) + 1]) / (2 * dy)));
                }
            }
            AInTime[n] = ATemp;
        }
        AInTimeGlobal = AInTime;

        // // write to file
        string filePath = Application.dataPath + @"\" +txtFileName+ ".txt";
        Debug.Log(filePath);
        using (StreamWriter writer = new StreamWriter(filePath))
        {
            for(int n = 0; n < AInTimeGlobal.Length; n++)
            {
                writer.WriteLine("");
                writer.WriteLine("TIME STEP N = " + n);
                for(int k=1;k<= vectToLin(xCount, yCount); k++)
                {
                    int[] lV = linToVect(k);
                    writer.WriteLine("k = "+k+" | i, j = "+lV[0]+", "+lV[1]+" | U, V = "+AInTimeGlobal[n].Values[2*(k-1)] + ", " + AInTimeGlobal[n].Values[2 * (k - 1)+1]);
                }
            }
        }
        for (int k = 1; k <= vectToLin(xCount, yCount); k++)
        {
            GameObject vV = vectorVisuals[k];
            vV.GetComponentInChildren<MeshRenderer>().material = new Material(white);
        }
    }
    void Update()
    {
        if (Time.time > 2)
        {
            int n = 0;
            n = ((int)((Time.time - 2) / dt*replaySpeed)) % nCount;
            Debug.Log(n*dt);

            for(int k = 1; k <= vectToLin(xCount, yCount); k++)
            {
                GameObject vV = vectorVisuals[k];

                float magnitude =Mathf.Max(-40,Mathf.Min(40,Mathf.Sqrt(Mathf.Pow((float)AInTimeGlobal[n].Values[2 * (k - 1)], 2) + Mathf.Pow((float)AInTimeGlobal[n].Values[2 * (k - 1) + 1], 2))));

                if (k == 50)
                {
                    Debug.Log(magnitude);
                }
                vV.transform.localScale = vectorScaleFactor*new Vector3(magnitude, 1, 1);
                vV.transform.rotation = Quaternion.Euler(0, 0, -xYToAngle(AInTimeGlobal[n].Values[3* (k - 1)], AInTimeGlobal[n].Values[3* (k - 1) + 1]));

                MeshRenderer[] mR = vV.GetComponentsInChildren<MeshRenderer>();
                for(int i=0; i < mR.Length; i++)
                {
                    mR[i].material.color= new Color(function_sigmoid(magnitude - 2), function_sigmoid(magnitude-2), function_sigmoid(magnitude - 2));
                }
            }
        }
    }
}