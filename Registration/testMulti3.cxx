#include <iostream>
#include <string>
#include <chrono>
#include "stdio.h"
using namespace std;

int main(int argc, char *argv[])
{
	int numberOfGridNodesInOneDimension[] = { 5 };
	int MaximumNumberOfIterations[] = { 20 };
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

			/*
			char* inputfolder = " C:/Users/sm5797/Dropbox/Research/ITK/ImageVolume/20160719/scan1/";
			char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/";
			char* labelMapOutputFolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan1/labelmap_Deformable/";
			char* labelMapInputFolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan1/labelmap_Deformable/";
			*/
			// Sanketh input folder
			char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/";
			char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROI/";
			char* labelMapOutputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/labelmap_Deformable/";
			char* labelMapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/";

			char outputfolder[300];
			sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/BSpline_GridNodeinOneDimension%d_Iterations%d_", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

			char fixedImageFile[300];
			sprintf_s(fixedImageFile, "%sminutes0.nii", fixedfolder);

			char movingImageFile[300];
			sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);
			//sprintf_s(movingImageFile, "%sSimilarity_IterationNum150_TransScale2_ScaleFactor5_Affine_IterationNum150RegistrationOutput.nii", inputfolder);
			char registeredImageFile[300];
			sprintf_s(registeredImageFile, "%sRegistrationOutput.nii", outputfolder);

			char differenceBeforeFile[300];
			sprintf_s(differenceBeforeFile, "%sdifferenceBeforeRegistration.nii.gz", outputfolder);

			char differenceAfterFile[300];
			sprintf_s(differenceAfterFile, "%sdifferenceAfterRegistration.nii.gz", outputfolder);

			char deformationField[300];
			sprintf_s(deformationField, "%sdeformationField.nrrd", outputfolder);

#pragma region labelMapFileName
			//output
			char labelmapOutput_SL[300];
			sprintf_s(labelmapOutput_SL, " %sSL_labelmap_registrationOutput.nii", outputfolder);

			char labelmapOutput_SM[300];
			sprintf_s(labelmapOutput_SM, " %sSM_labelmap_registrationOutput.nii", outputfolder);
			char labelmapOutput_SM_e1[300];
			sprintf_s(labelmapOutput_SM_e1, " %sSM_e1_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_SM_e2[300];
			sprintf_s(labelmapOutput_SM_e2, " %sSM_e2_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_SM_e3[300];
			sprintf_s(labelmapOutput_SM_e3, " %sSM_e3_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_SM_e4[300];
			sprintf_s(labelmapOutput_SM_e4, " %sSM_e4_labelmap_registrationOutput.nii", labelMapOutputFolder);

			char labelmapOutput_ST[300];
			sprintf_s(labelmapOutput_ST, " %sST_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_ST_e1[300];
			sprintf_s(labelmapOutput_ST_e1, " %sST_e1_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_ST_e2[300];
			sprintf_s(labelmapOutput_ST_e2, " %sST_e2_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_ST_e3[300];
			sprintf_s(labelmapOutput_ST_e3, " %sST_e3_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_ST_e4[300];
			sprintf_s(labelmapOutput_ST_e4, " %sST_e4_labelmap_registrationOutput.nii", labelMapOutputFolder);

			char labelmapOutput_SV[300];
			sprintf_s(labelmapOutput_SV, " %sSV_labelmap_registrationOutput.nii", outputfolder);
			char labelmapOutput_SV_e1[300];
			sprintf_s(labelmapOutput_SV_e1, " %sSV_e1_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_SV_e2[300];
			sprintf_s(labelmapOutput_SV_e2, " %sSV_e2_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_SV_e3[300];
			sprintf_s(labelmapOutput_SV_e3, " %sSV_e3_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_SV_e4[300];
			sprintf_s(labelmapOutput_SV_e4, " %sSV_e4_labelmap_registrationOutput.nii", labelMapOutputFolder);
			char labelmapOutput_SV_e5[300];
			sprintf_s(labelmapOutput_SV_e5, " %sSV_e5_labelmap_registrationOutput.nii", labelMapOutputFolder);

			//Input
			char labelmapInput_SV[300];
			sprintf_s(labelmapInput_SV, " %satlasVolume_SV_labelmap.nii", labelMapInputFolder);
			char labelmapInput_SV_e1[300];
			sprintf_s(labelmapInput_SV_e1, " %satlasVolume_SV_labelmap_e1.nii", labelMapInputFolder);
			char labelmapInput_SV_e2[300];
			sprintf_s(labelmapInput_SV_e2, " %satlasVolume_SV_labelmap_e2.nii", labelMapInputFolder);
			char labelmapInput_SV_e3[300];
			sprintf_s(labelmapInput_SV_e3, " %satlasVolume_SV_labelmap_e3.nii", labelMapInputFolder);
			char labelmapInput_SV_e4[300];
			sprintf_s(labelmapInput_SV_e4, " %satlasVolume_SV_labelmap_e4.nii", labelMapInputFolder);
			char labelmapInput_SV_e5[300];
			sprintf_s(labelmapInput_SV_e5, " %satlasVolume_SV_labelmap_e5.nii", labelMapInputFolder);

			char labelmapInput_ST[300];
			sprintf_s(labelmapInput_ST, " %satlasVolume_ST_labelmap.nii", labelMapInputFolder);
			char labelmapInput_ST_e1[300];
			sprintf_s(labelmapInput_ST_e1, " %satlasVolume_ST_labelmap_e1.nii", labelMapInputFolder);
			char labelmapInput_ST_e2[300];
			sprintf_s(labelmapInput_ST_e2, " %satlasVolume_ST_labelmap_e2.nii", labelMapInputFolder);
			char labelmapInput_ST_e3[300];
			sprintf_s(labelmapInput_ST_e3, " %satlasVolume_ST_labelmap_e3.nii", labelMapInputFolder);
			char labelmapInput_ST_e4[300];
			sprintf_s(labelmapInput_ST_e4, " %satlasVolume_ST_labelmap_e4.nii", labelMapInputFolder);

			char labelmapInput_SM[300];
			sprintf_s(labelmapInput_SM, " %satlasVolume_SM_labelmap.nii", labelMapInputFolder);
			char labelmapInput_SM_e1[300];
			sprintf_s(labelmapInput_SM_e1, " %satlasVolume_SM_labelmap_e1.nii", labelMapInputFolder);
			char labelmapInput_SM_e2[300];
			sprintf_s(labelmapInput_SM_e2, " %satlasVolume_SM_labelmap_e2.nii", labelMapInputFolder);
			char labelmapInput_SM_e3[300];
			sprintf_s(labelmapInput_SM_e3, " %satlasVolume_SM_labelmap_e3.nii", labelMapInputFolder);
			char labelmapInput_SM_e4[300];
			sprintf_s(labelmapInput_SM_e4, " %satlasVolume_SM_labelmap_e4.nii", labelMapInputFolder);

			char labelmapInput_SL[300];
			sprintf_s(labelmapInput_SL, " %sSL_labelmap_registrationOutput.nii", labelMapInputFolder);
#pragma endregion

			char parameter[50];
			sprintf_s(parameter, " %d %d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

			char command[8000] = "";
			strcat_s(command, "E:/Sanketh_Moudgalya_Directory/Nate/C++_codes/built_codes/Regi/Debug/MultiStage3.exe ");

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

			strcat_s(command, labelmapInput_SV_e1);
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
			strcat_s(command, labelmapOutput_SM_e4);

			cout << command << endl;
			std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

			system(command);
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

			std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes" << std::endl;

		}
	}

	//system("pause");
	return EXIT_SUCCESS;

}