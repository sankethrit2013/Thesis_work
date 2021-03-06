#include <iostream>
#include <string>
#include "stdio.h"
#include <chrono>

using namespace std;

int main()
{
	double translationFactor[] = { 0.1};
	int iterationsNumber[] = { 300 };
	double learningRate[] = { 0.1 };
	double minimumStepLength[] = { 0.00001 };
	for (int t = 2; t <= 4; t++)
	{
		if (t != 1)
		{
			//sanketh code
			std::cout << "scan" << t << " computing" << std::endl;

			for (int l = 0; l < sizeof(minimumStepLength) / sizeof(minimumStepLength[0]); l++)
			{
				for (int j = 0; j < sizeof(learningRate) / sizeof(learningRate[0]); j++)
				{
					for (int k = 0; k < sizeof(iterationsNumber) / sizeof(iterationsNumber[0]); k++)
					{
						for (int i = 0; i < sizeof(translationFactor) / sizeof(translationFactor[0]); i++)
						{
							//char* scansfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/";

							//char* labelmapInputFolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan1/labelmap_Deformable_Grid10_Iter100/";
							//char labelmapOutputFolder[300];
							//sprintf_s(labelmapOutputFolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan%d/labelmap_Deformable_Grid10_Iter100_fromscan1/", t);

							//char fixedImageFile[300];
							//sprintf_s(fixedImageFile, "%sscan%d_leftEar.nii", scansfolder, t);
							//char movingImageFile[300];
							//sprintf_s(movingImageFile, "%sscan1_leftEar.nii", scansfolder, t);

							//char ResampleInputFile[300];
							///*sprintf_s(ResampleInputFile, "%sscan1/BSpline_GridNodeinOneDimension8_Iterations50_RegistrationOutput.nii", scansfolder);*/
							//sprintf_s(ResampleInputFile, "%sscan1/BSpline_GridNodeinOneDimension10_Iterations100_RegistrationOutput.nii", scansfolder);
							//char ResampleOutputFile[300];
							///*sprintf_s(ResampleOutputFile, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension8_Iterations50_RegistrationOutput.nii", scansfolder, t);*/
							//sprintf_s(ResampleOutputFile, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension10_Iterations100_RegistrationOutput.nii", scansfolder, t);


							//char differenceAfterFile_scan2scan[300];
							//sprintf_s(differenceAfterFile_scan2scan, "%sscan%d/Rigid_toScan1_iterationNum%d_translationFactor%d_learningRate%g_minStepLength%g_differenceAfterRegistration.nii", scansfolder, t, iterationsNumber[k], translationFactor[i], learningRate[j], minimumStepLength[l]);

							//char differenceAfterFile_scan2atlas[300];
							///*sprintf_s(differenceAfterFile_scan2atlas, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension8_Iterations50_differenceAfterRegistration.nii", scansfolder, t);*/
							//sprintf_s(differenceAfterFile_scan2atlas, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension10_Iterations100_differenceAfterRegistration.nii", scansfolder, t);


							//Sanketh folders
							//20160906 scans
							/*
							char* scansfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/";

							char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/scan1/labelmap_Deformable_Grid10_Iter100/labelmap_Deformable_Grid10_Iter100";
							char labelmapOutputFolder[300];
							sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/scan%d/labelmap_Deformable_Grid10_Iter100_fromscan1/", t);
							*/

							/*-----------------------------------------------------------------------------------------------------------------*/
							//20160629 scans
							
							/*char* scansfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/";

							char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan1/labelmap_Deformable_Grid10_Iter100/";
							char labelmapOutputFolder[300];
							sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan%d/labelmap_Deformable_Grid10_Iter100_fromscan1/", t);*/
							/*-----------------------------------------------------------------------------------------------------------------*/

							///*-----------------------------------------------------------------------------------------------------------------*/
							////20170111 scans

							//char* scansfolder = " e:/sanketh_moudgalya_directory/nate/imagevolume/20170123/";
							//char* baselinefolder = " e:/sanketh_moudgalya_directory/nate/imagevolume/20170123/scan1/";
							//char otherscansfolder[300];
							//sprintf_s(otherscansfolder, " e:/sanketh_moudgalya_directory/nate/imagevolume/20170123/scan%d/", t);


							//char* labelmapinputfolder = " e:/sanketh_moudgalya_directory/nate/imagevolume/20170123/scan1/labelmap_deformable_grid10_iter150/";
							//char labelmapoutputfolder[300];
							//sprintf_s(labelmapoutputfolder, " e:/sanketh_moudgalya_directory/nate/imagevolume/20170123/scan%d/labelmap_deformable_grid10_iter150_fromscan1/", t);
							///*-----------------------------------------------------------------------------------------------------------------*/

							///*-----------------------------------------------------------------------------------------------------------------*/
							////20170606 scans

							//char* scansfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/";
							//char* baselinefolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/";
							//char otherScansFolder[300];
							//sprintf_s(otherScansFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan%d/", t);


							//char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/finalLM/";
							//char labelmapOutputFolder[300];
							//sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan%d/labelmap_deformable_grid10_iter100_fromscan1/", t);
							///*-----------------------------------------------------------------------------------------------------------------*/

							///*-----------------------------------------------------------------------------------------------------------------*/
							////20170619 scans

							//char* scansfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/";
							//char* baselinefolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan1/";
							//char otherScansFolder[300];
							//sprintf_s(otherScansFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan%d/", t);


							//char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/BSpline/";
							//char labelmapOutputFolder[300];
							//sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan%d/labelmap_deformable_grid10_iter100_fromscan1/", t);
							///*-----------------------------------------------------------------------------------------------------------------*/

							/*-----------------------------------------------------------------------------------------------------------------*/
							//20171017 scans

							char* scansfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/";
							char* baselinefolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan1/";
							char otherScansFolder[300];
							sprintf_s(otherScansFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan%d/", t);


							char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/BSpline/";
							char labelmapOutputFolder[300];
							sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan%d/labelmap_deformable_grid10_iter300_fromscan1/", t);
							/*-----------------------------------------------------------------------------------------------------------------*/

							char fixedImageFile[300];
							sprintf_s(fixedImageFile, "%sscan%d.nii", otherScansFolder, t);
							char movingImageFile[300];
							sprintf_s(movingImageFile, "%sscan1.nii", baselinefolder);

							char ResampleInputFile[300];
							sprintf_s(ResampleInputFile, "%sscan1/BSpline_GridNodeinOneDimension10_Iterations300RegistrationOutput.nii", scansfolder);
							//sprintf_s(ResampleInputFile, "%sscan1/BSpline_GridNodeinOneDimension10_Iterations100RegistrationOutput.nii", scansfolder);
							char ResampleOutputFile[300];
							sprintf_s(ResampleOutputFile, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension10_Iterations300RegistrationOutput.nii", scansfolder, t);
							//sprintf_s(ResampleOutputFile, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension10_Iterations100_RegistrationOutput.nii", scansfolder, t);
							//sprintf_s(ResampleOutputFile, "%sRegistrationOutput.nii", scansfolder);


							char differenceAfterFile_scan2scan[300];
							sprintf_s(differenceAfterFile_scan2scan, "%sscan%d/Rigid_toScan1_iterationNum%d_translationFactor%g_learningRate%g_minStepLength%g_differenceAfterRegistration.nii", scansfolder, t, iterationsNumber[k], translationFactor[i], learningRate[j], minimumStepLength[l]);

							char differenceAfterFile_scan2atlas[300];
							/*sprintf_s(differenceAfterFile_scan2atlas, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension8_Iterations50_differenceAfterRegistration.nii", scansfolder, t);*/
							sprintf_s(differenceAfterFile_scan2atlas, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension10_Iterations300_differenceAfterRegistration.nii", scansfolder, t);

							char parameter[50];
							sprintf_s(parameter, " %g %d %g %g", translationFactor[i], iterationsNumber[k], learningRate[j], minimumStepLength[l]);
#pragma region labelMapFileName
							//label map output
							char labelmapOutput_SL[300];
							sprintf_s(labelmapOutput_SL, " %sSL_labelmap_Iter%d_registrationOutput.nii", labelmapOutputFolder, iterationsNumber[k]);

							char labelmapOutput_SM[300];
							sprintf_s(labelmapOutput_SM, " %sSM_labelmap_Iter%d_registrationOutput.nii", labelmapOutputFolder, iterationsNumber[k]);

							char labelmapOutput_ST[300];
							sprintf_s(labelmapOutput_ST, " %sST_labelmap_Iter%d_registrationOutput.nii", labelmapOutputFolder, iterationsNumber[k]);
							
							char labelmapOutput_SV[300];
							sprintf_s(labelmapOutput_SV, " %sSV_labelmap_Iter%d_registrationOutput.nii", labelmapOutputFolder, iterationsNumber[k]);
							
							//label map input
							char labelmapInput_SL[300];
							sprintf_s(labelmapInput_SL, " %slabelmap_Deformable_Grid10_Iter300SL_labelmap_registrationOutput.nii", labelmapInputFolder);

							char labelmapInput_SM[300];
							sprintf_s(labelmapInput_SM, " %slabelmap_Deformable_Grid10_Iter300SM_labelmap_registrationOutput.nii", labelmapInputFolder);
							
							char labelmapInput_ST[300];
							sprintf_s(labelmapInput_ST, " %slabelmap_Deformable_Grid10_Iter300ST_labelmap_registrationOutput.nii", labelmapInputFolder);
							
							char labelmapInput_SV[300];
							sprintf_s(labelmapInput_SV, " %slabelmap_Deformable_Grid10_Iter300SV_labelmap_registrationOutput.nii", labelmapInputFolder);
							
#pragma endregion
							char command[8000] = "";

							strcat_s(command, "E:/Sanketh_Moudgalya_Directory/Nate/C++_codes/built_codes/Regi/Debug/RegistrationRigid3D.exe ");
							strcat_s(command, fixedImageFile);
							strcat_s(command, movingImageFile);
							strcat_s(command, ResampleInputFile);

							strcat_s(command, ResampleOutputFile);
							strcat_s(command, differenceAfterFile_scan2scan);
							strcat_s(command, differenceAfterFile_scan2atlas);
							strcat_s(command, parameter);
							strcat_s(command, labelmapInput_SV); //argv[11]
							strcat_s(command, labelmapOutput_SV);
							strcat_s(command, labelmapInput_ST);
							strcat_s(command, labelmapOutput_ST);
							strcat_s(command, labelmapInput_SM);
							strcat_s(command, labelmapOutput_SM);
							strcat_s(command, labelmapInput_SL);
							strcat_s(command, labelmapOutput_SL);

							cout << command << endl;
							std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

							system(command);
							std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

							std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes" << std::endl;

						}
					}
				}
			}
		}
	}
	return 0;
}