#include <iostream>
#include <string>
#include <chrono>
#include "stdio.h"
using namespace std;

int main(int argc, char *argv[])
{
	//int numberOfGridNodesInOneDimension[] = {8, 10};
	int numberOfGridNodesInOneDimension[] = {10};
	//int MaximumNumberOfIterations[] = {50, 100};
	int MaximumNumberOfIterations[] = { 100 };
	for (size_t t = 0; t < 5; t++)
	{
		if (t == 1)
		{
			for (int j = 0; j < sizeof(MaximumNumberOfIterations) / sizeof(MaximumNumberOfIterations[0]); j++)
			{
				for (int k = 0; k < sizeof(numberOfGridNodesInOneDimension) / sizeof(numberOfGridNodesInOneDimension[0]); k++)
				{
					//setting for reg using prevvious stage output
					//char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20150801/Affine/Scale1_0.2_IterationNum500_LearningRate1_";
					//char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20150801";
					//char outputfolder[300];
					//sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20150801/BSpline/GridNodeinOneDimension%d_Iterations%d_",numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					//char fixedImageFile[300];
					//sprintf_s(fixedImageFile, "%s/dicom_ear_male_leftEar.nii", fixedfolder);
					//char movingImageFile[300];
					//sprintf_s(movingImageFile, "%sRegistrationOutput", inputfolder);

					/*char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/Affine_0108first_atlasnew3";
					char* outputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/Affine_0108first_atlasnew3";
					char fixedImageFile[300];
					sprintf_s(fixedImageFile, "%s/first_base_line_ear.nii.gz", inputfolder);*/

					//Zhenlin inputs
					/*
					char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/";
					char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/";
					char labelmapOutputFolder[300];
					sprintf_s(labelmapOutputFolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan%d/labelmap_Deformable_Grid%d_Iter%d/", t, numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					char outputfolder[300];
					sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan%d/BSpline_GridNodeinOneDimension%d_Iterations%d_", t, numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);
					*/

					/*-------------------------------------------------------------------------------------------------*/
					/*20160906 scans*/
					//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROI/";
					//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/MultiStage/";
					////// Sanketh addition
					//char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/";

					//char labelmapOutputFolder[300];
					//////sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/Bspline/labelmap_Deformable_Grid%d_Iter%d/", t, numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);
					//sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/Bspline/labelmap_Deformable_Grid%d_Iter%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					//char outputfolder[300];
					//////sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan%d/BSpline_GridNodeinOneDimension%d_Iterations%d_", t, numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);
					//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/Bspline/BSpline_GridNodeinOneDimension%d_Iterations%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					/*------------------------------------------------------------------------------------------------*/
					// 20160629 scans

					//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/";
					//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/MultiStage/";
					//// Sanketh addition
					//char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/";

					//char labelmapOutputFolder[300];
					////sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/Bspline/labelmap_Deformable_Grid%d_Iter%d/", t, numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);
					//sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/Bspline/labelmap_Deformable_Grid%d_Iter%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					//char outputfolder[300];
					////sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan%d/BSpline_GridNodeinOneDimension%d_Iterations%d_", t, numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);
					//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan1/BSpline_GridNodeinOneDimension%d_Iterations%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);
					/*------------------------------------------------------------------------------------------------*/

					/*-----------------------------------------------------------------------------------------------------*/
					///*20170111 scans*/
					//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/";
					//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/scan1/";
	
					//char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/";

					//char labelmapOutputFolder[300];
					//sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/Bspline/labelmap_Deformable_Grid%d_Iter%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					//char outputfolder[300];
					//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/scan1/BSpline_GridNodeinOneDimension%d_Iterations%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					///*-----------------------------------------------------------------------------------------------------*/

					///*20170123 scans*/
					//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/";
					//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/";

					//char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/";

					//char labelmapOutputFolder[300];
					//sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/Bspline/labelmap_Deformable_Grid%d_Iter%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					//char outputfolder[300];
					//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/BSpline_GridNodeinOneDimension%d_Iterations%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					///*-----------------------------------------------------------------------------------------------------*/

					///*-----------------------------------------------------------------------------------------------------*/

					///*20170606 scans*/
					//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/";
					//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/";

					//char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/";

					//char labelmapOutputFolder[300];
					//sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/Bspline/labelmap_Deformable_Grid%d_Iter%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					//char outputfolder[300];
					//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/BSpline_GridNodeinOneDimension%d_Iterations%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					///*-----------------------------------------------------------------------------------------------------*/

					///*-----------------------------------------------------------------------------------------------------*/

					///*20170612 scans*/
					//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/";
					//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan1/";

					//char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/";

					//char labelmapOutputFolder[300];
					//sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/Bspline/labelmap_Deformable_Grid%d_Iter%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					//char outputfolder[300];
					//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan1/BSpline_GridNodeinOneDimension%d_Iterations%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					///*-----------------------------------------------------------------------------------------------------*/

					/*-----------------------------------------------------------------------------------------------------*/

					/*20170606 scans*/
					char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/";
					char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/";

					char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/";

					char labelmapOutputFolder[300];
					sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/FinalLM/labelmap_Deformable_Grid%d_Iter%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					char outputfolder[300];
					sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/BSpline_GridNodeinOneDimension%d_Iterations%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					/*-----------------------------------------------------------------------------------------------------*/


					char fixedImageFile[300];
					sprintf_s(fixedImageFile, "%sscan1_other.nii", fixedfolder);
					char movingImageFile[300];
					sprintf_s(movingImageFile, "%satlas_isotropic01638.nii", inputfolder);
					//sprintf_s(movingImageFile, "%sbaseline_Similarity_IterationNum150_TransScale2_ScaleFactor5_Affine_IterationNum100RegistrationOutput.nii", fixedfolder);

					char registeredImageFile[300];
					sprintf_s(registeredImageFile, "%sRegistrationOutput.nii", outputfolder);

					char differenceBeforeFile[300];
					sprintf_s(differenceBeforeFile, "%sdifferenceBeforeRegistration.nii.gz", outputfolder);

					char differenceAfterFile[300];
					sprintf_s(differenceAfterFile, "%sdifferenceAfterRegistration.nii.gz", outputfolder);

					char deformationField[300];
					sprintf_s(deformationField, "%sdeformationField.nrrd", outputfolder);
					/*-----------------------------------------------------------------------------------------------------*/

#pragma region labelMapFileName
					//output
					char labelmapOutput_SL[200];
					sprintf_s(labelmapOutput_SL, " %sSL_labelmap_registrationOutput.nii", labelmapOutputFolder);

					char labelmapOutput_SM[200];
					sprintf_s(labelmapOutput_SM, " %sSM_labelmap_registrationOutput.nii", labelmapOutputFolder);
					/*char labelmapOutput_SM_e1[200];
					sprintf_s(labelmapOutput_SM_e1, " %sSM_e1_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_SM_e2[200];
					sprintf_s(labelmapOutput_SM_e2, " %sSM_e2_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_SM_e3[200];
					sprintf_s(labelmapOutput_SM_e3, " %sSM_e3_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_SM_e4[200];
					sprintf_s(labelmapOutput_SM_e4, " %sSM_e4_labelmap_registrationOutput.nii", labelmapOutputFolder);*/

					char labelmapOutput_ST[200];
					sprintf_s(labelmapOutput_ST, " %sST_labelmap_registrationOutput.nii", labelmapOutputFolder);
					/*char labelmapOutput_ST_e1[200];
					sprintf_s(labelmapOutput_ST_e1, " %sST_e1_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_ST_e2[200];
					sprintf_s(labelmapOutput_ST_e2, " %sST_e2_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_ST_e3[200];
					sprintf_s(labelmapOutput_ST_e3, " %sST_e3_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_ST_e4[200];
					sprintf_s(labelmapOutput_ST_e4, " %sST_e4_labelmap_registrationOutput.nii", labelmapOutputFolder);*/

					char labelmapOutput_SV[200];
					sprintf_s(labelmapOutput_SV, " %sSV_labelmap_registrationOutput.nii", labelmapOutputFolder);
					/*char labelmapOutput_SV_e1[200];
					sprintf_s(labelmapOutput_SV_e1, " %sSV_e1_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_SV_e2[200];
					sprintf_s(labelmapOutput_SV_e2, " %sSV_e2_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_SV_e3[200];
					sprintf_s(labelmapOutput_SV_e3, " %sSV_e3_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_SV_e4[200];
					sprintf_s(labelmapOutput_SV_e4, " %sSV_e4_labelmap_registrationOutput.nii", labelmapOutputFolder);
					char labelmapOutput_SV_e5[200];
					sprintf_s(labelmapOutput_SV_e5, " %sSV_e5_labelmap_registrationOutput.nii", labelmapOutputFolder);*/
					char labelmapOutput_CA[300];
					sprintf_s(labelmapOutput_CA, " %sCA_labelmap_registrationOutput.nii", labelmapOutputFolder);
					
					/*char labelmapOutput_Comb[300];
					sprintf_s(labelmapOutput_Comb, " %slabelmap_combined01638_registrationOutput.nii", labelmapOutputFolder);*/


					//Input
					//Sanketh - change input folder to labelmapInputFolder
					char labelmapInput_SV[200];
					sprintf_s(labelmapInput_SV, " %satlasVolume_SV_labelmap.nii", labelmapInputFolder);
					/*char labelmapInput_SV_e1[200];
					sprintf_s(labelmapInput_SV_e1, " %satlasVolume_SV_labelmap_e1.nii", labelmapInputFolder);
					char labelmapInput_SV_e2[200];
					sprintf_s(labelmapInput_SV_e2, " %satlasVolume_SV_labelmap_e2.nii", labelmapInputFolder);
					char labelmapInput_SV_e3[200];
					sprintf_s(labelmapInput_SV_e3, " %satlasVolume_SV_labelmap_e3.nii", labelmapInputFolder);
					char labelmapInput_SV_e4[200];
					sprintf_s(labelmapInput_SV_e4, " %satlasVolume_SV_labelmap_e4.nii", labelmapInputFolder);
					char labelmapInput_SV_e5[200];
					sprintf_s(labelmapInput_SV_e5, " %satlasVolume_SV_labelmap_e5.nii", labelmapInputFolder);*/

					char labelmapInput_ST[200];
					sprintf_s(labelmapInput_ST, " %satlasVolume_ST_labelmap.nii", labelmapInputFolder);
					/*char labelmapInput_ST_e1[200];
					sprintf_s(labelmapInput_ST_e1, " %satlasVolume_ST_labelmap_e1.nii", labelmapInputFolder);
					char labelmapInput_ST_e2[200];
					sprintf_s(labelmapInput_ST_e2, " %satlasVolume_ST_labelmap_e2.nii", labelmapInputFolder);
					char labelmapInput_ST_e3[200];
					sprintf_s(labelmapInput_ST_e3, " %satlasVolume_ST_labelmap_e3.nii", labelmapInputFolder);
					char labelmapInput_ST_e4[200];
					sprintf_s(labelmapInput_ST_e4, " %satlasVolume_ST_labelmap_e4.nii", labelmapInputFolder);*/

					char labelmapInput_SM[200];
					sprintf_s(labelmapInput_SM, " %satlasVolume_SM_labelmap.nii", labelmapInputFolder);
					/*char labelmapInput_SM_e1[200];
					sprintf_s(labelmapInput_SM_e1, " %satlasVolume_SM_labelmap_e1.nii", labelmapInputFolder);
					char labelmapInput_SM_e2[200];
					sprintf_s(labelmapInput_SM_e2, " %satlasVolume_SM_labelmap_e2.nii", labelmapInputFolder);
					char labelmapInput_SM_e3[200];
					sprintf_s(labelmapInput_SM_e3, " %satlasVolume_SM_labelmap_e3.nii", labelmapInputFolder);
					char labelmapInput_SM_e4[200];
					sprintf_s(labelmapInput_SM_e4, " %satlasVolume_SM_labelmap_e4.nii", labelmapInputFolder);*/

					char labelmapInput_SL[200];
					sprintf_s(labelmapInput_SL, " %satlasVolume_SL_labelmap.nii", labelmapInputFolder);
					char labelmapInput_CA[300];
					sprintf_s(labelmapInput_CA, " %sCA_labelmap.nii", labelmapInputFolder);
					
					/*char labelmapInput_Comb[300];
					sprintf_s(labelmapInput_Comb, " %slabelmap_combined01638.nii", labelmapInputFolder);*/
#pragma endregion

					char parameter[50];
					sprintf_s(parameter, " %d %d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					char command[8000] = "";
					strcat_s(command, "E:/Sanketh_Moudgalya_Directory/Nate/C++_codes/built_codes/Regi/Debug/DeformableRegistrationBSpline3D.exe ");

					strcat_s(command, fixedImageFile);
					strcat_s(command, movingImageFile);
					strcat_s(command, registeredImageFile);
					strcat_s(command, differenceAfterFile);
					strcat_s(command, differenceBeforeFile);
					strcat_s(command, parameter);
					strcat_s(command, deformationField);
					strcat_s(command, labelmapInput_SV);
					strcat_s(command, labelmapOutput_SV);
					strcat_s(command, labelmapInput_ST);
					strcat_s(command, labelmapOutput_ST);
					strcat_s(command, labelmapInput_SM);
					strcat_s(command, labelmapOutput_SM);
					strcat_s(command, labelmapInput_SL);
					strcat_s(command, labelmapOutput_SL);
					strcat_s(command, labelmapInput_CA);
					strcat_s(command, labelmapOutput_CA);


					/*strcat_s(command, labelmapInput_SV_e1);
					strcat_s(command, labelmapOutput_SV_e1);
					strcat_s(command, labelmapInput_SV_e2);
					strcat_s(command, labelmapOutput_SV_e2);
					strcat_s(command, labelmapInput_SV_e3);
					strcat_s(command, labelmapOutput_SV_e3);
					strcat_s(command, labelmapInput_SV_e4);
					strcat_s(command, labelmapOutput_SV_e4);
					strcat_s(command, labelmapInput_SV_e5);
					strcat_s(command, labelmapOutput_SV_e5);

					strcat_s(command, labelmapInput_ST_e1);
					strcat_s(command, labelmapOutput_ST_e1);
					strcat_s(command, labelmapInput_ST_e2);
					strcat_s(command, labelmapOutput_ST_e2);
					strcat_s(command, labelmapInput_ST_e3);
					strcat_s(command, labelmapOutput_ST_e3);
					strcat_s(command, labelmapInput_ST_e4);
					strcat_s(command, labelmapOutput_ST_e4);

					strcat_s(command, labelmapInput_SM_e1);
					strcat_s(command, labelmapOutput_SM_e1);
					strcat_s(command, labelmapInput_SM_e2);
					strcat_s(command, labelmapOutput_SM_e2);
					strcat_s(command, labelmapInput_SM_e3);
					strcat_s(command, labelmapOutput_SM_e3);
					strcat_s(command, labelmapInput_SM_e4);
					strcat_s(command, labelmapOutput_SM_e4);*/

					cout << command << endl;
					std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

					system(command);
					std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

					std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes" << std::endl;

				}
			}

		}
	}
	return EXIT_SUCCESS;
}