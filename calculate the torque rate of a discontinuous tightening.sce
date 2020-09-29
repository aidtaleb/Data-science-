//============================================================================
// nom : caracterisation d'assemblage en vissage discontinu.sce
// auteur : ID TALEB ABDERRAHIM 
// date de création : 2017-05-15
// dates de modification : ...
//----------------------------------------------------------------------------
// version de Scilab : 5.4.1
// module Atoms requis : aucun
//----------------------------------------------------------------------------
// Objectif : Déterminer la nature de l'assemblage.
// Entrées : courbes de vissage 
// Sorties : torque rate, modélisation de l'assemblage
//------------------------------------------------------------------------------
clear();
clc(); 
//=============================================================================
// importation et reception des fichiers texts 
//=============================================================================

stacksize('max')

global Stop_programme;
global Run_Analyze;
global TableMessage;
global Actualiser;
//==============================================================================
//                         LES FONCTIONS 
//==============================================================================
//========================================================================//
// this function "Torquerate_method_median" is used to calculate the torque rate in continuous domaine//
//and also discontinuous domaine after the estimation of the points       //
//========================================================================//
//==============================================================================
function [Torque_rate]=Torquerate_method_median(Angle,Torque,Angle_span,Taille_du_Buffer)

    s=0;
    index=0;
    j=1;
    median_Angle=0;
    median_Torque=0;
    indice_last_angle_step=1;// ns permet de conserver l'indice du i où moment on a arrivé à l'angle step
    Increment_buffer=1;
    Indice_first_buffer=1;
    Increment_Torque_rate=1;
    for i=1:length(Angle ) 
        // s'il s'agit du premiere indice on le stoke dans notre buffer 
        if i==1 then 
            Buffer_angle(Increment_buffer)=Angle (i)-Angle(1);
            Buffer_torque(Increment_buffer)=Torque(i);
            Indice_end_buffer =Increment_buffer;
            Increment_buffer=Increment_buffer+1;
            // sinon on paase au suivant 
        else 
            // on teste que la premiere difference entre le premier élement et celui de l'étape i est supérieur ou égale à 
            //l'angle step si ou onstoke l'element dans le buffer et on garde l'indice de l'élemnt dans la memoire

            if Angle(i)- Angle (indice_last_angle_step)>= Angle_span /(Taille_du_Buffer) // est que la difference est superieur ou égale à l'angle step
                Buffer_angle(Increment_buffer)=Angle (i)-Angle (1);                  //stokage de'angle dans le buffer 
                Buffer_torque(Increment_buffer)=Torque (i); 
                // tabtorqurate(Increment_buffer)  =Torque_rate_measured(i);                                //stokage du couple correspondant 
                indice_last_angle_step =i;                                                                // stokage de l'indice de l'élement qui correspond à l'angle step
                Indice_end_buffer=Increment_buffer;                                                               // stokage d'indice de dernier élement stoké dans le tableau
                Increment_buffer=Increment_buffer+1;
            end 

            if Buffer_angle(Indice_end_buffer)-Buffer_angle(Indice_first_buffer)>= Angle_span then                         // est ce que la difference est superieur ou égale à l'angle span
                //boucle pour chercher le melieu du tableau
                pp=Indice_first_buffer;                                                              //on stok l'indice du premier élement du buffer en memoire 
                while s< Angle_span/2
                    s=s+ Buffer_angle(pp+1)-Buffer_angle(pp);
                    Indice_median_buffer=pp;
                    pp=pp+1;
                end 
                s=0; // remise de lindice s à 0
                // calcule du torque rate 
                DeltAngle= mean(Buffer_angle(Indice_first_buffer:Indice_median_buffer))- mean(Buffer_angle(Indice_median_buffer+1:Indice_end_buffer));
                DeltTorque= mean(Buffer_torque(Indice_first_buffer:Indice_median_buffer))- mean(Buffer_torque(Indice_median_buffer+1:Indice_end_buffer));
                if DeltAngle <> 0 then 
                    Torque_rate(Increment_Torque_rate)=  DeltTorque/ DeltAngle;

                    Increment_Torque_rate=Increment_Torque_rate+1;                                                                 // indice d'incrementation du tableau torque rate 
                else 
                    Torque_rate(Increment_Torque_rate)= 0;
                    Increment_Torque_rate=Increment_Torque_rate+1; 
                end 
                Indice_first_buffer = Indice_first_buffer+1;             // déplacement au deuxieme élement du buffer      
            end 
        end


    end

endfunction

//==============================================================================
//this function "Torque_rate_discontinu" calculates the torque rate en vissage discontinue 
//==============================================================================
function [Torque_rate]=Torque_rate_discontinu(Angle,Torque)
    j=1;
    for i=1:length(Angle) 
        if i>1 then 
            if Angle(i)-Angle(i-1)<>0 then 
                Torque_rate(j)=(Torque(i)-Torque(i-1))/(Angle(i)-Angle(i-1));
                j=j+1;
            else 
                if j>1 then 
                    Torque_rate(j)= Torque_rate(j-1);
                    j=j+1;
                else
                    Torque_rate(j)=0;
                    j=j+1;
                end

            end 
        end
    end

endfunction
//==============================================================================
// this faunction allows us to know the assembly is multi-stiffness or not 
//==============================================================================
function [ count, median_angle_1,median_torque_1, median_angle_2,median_torque_2]=poly_stiffness(Angle,Torque,Torque_rate)

    for i=1: length(Torque_rate)
        if i>1 then 
            if abs(Torque_rate(i)- Torque_rate(i-1))>0.7 then 
                median_angle_1=Angle(i-1);
                median_angle_2=Angle(i+1);
                median_torque_1=Torque(i-1);
                median_torque_2=Torque(i+1);
                count=1; 
                break;
            else 
                median_angle_1=0;
                median_angle_2=0;
                median_torque_1=0;
                median_torque_2=0;
                count=0;
            end
        end
    end

endfunction
//==============================================================================
//        initialization of variable 
//==============================================================================
Stop_programme = 0;
Run_Analyze = 0;
TableMessage = ["  >","  >","  >","  >","  >","  >","  >","  >","  >"];
Actualiser=0;
//==============================================================================
//this function "abort_callback" closes the programme
//==============================================================================
function abort_callback(handles)
    global Stop_programme;

    Stop_programme = 1;
    close(f);
    close(w);

    //clear
    abort
endfunction
//==============================================================================
// this function "test_callback" starts the programme 
//==============================================================================
function test_callback(handles)
    global Run_Analyze;
    Run_Analyze = 1;
endfunction
//==============================================================================
//this function "Actualisation_callback" update the programme
//==============================================================================
function Actualisation_callback (handles)
    global Actualiser;
    Actualiser=1;

endfunction
//==============================================================================
//this function "Message_console" creates the message console
//==============================================================================
function Message_console(handles,MyMessage)
    global TableMessage;
    TableMessage(2:9) = TableMessage(1:8);
    TableMessage(1) = msprintf('  > %s',MyMessage);
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(9),'Position',[10 460 800 20],'fontsize',14,'FontWeight','demi');
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(8),'Position',[10 440 800 20],'fontsize',14,'FontWeight','demi');
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(7),'Position',[10 420 800 20],'fontsize',14,'FontWeight','demi');
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(6),'Position',[10 400 800 20],'fontsize',14,'FontWeight','demi');
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(5),'Position',[10 380 800 20],'fontsize',14,'FontWeight','demi');
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(4),'Position',[10 360 800 20],'fontsize',14,'FontWeight','demi');
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(3),'Position',[10 340 800 20],'fontsize',14,'FontWeight','demi');
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(2),'Position',[10 320 800 20],'fontsize',14,'FontWeight','demi');
    handles.message_consol=uicontrol(w,"style","text","string",TableMessage(1),'Position',[10 300 800 20],'fontsize',14,'FontWeight','demi');
endfunction
//==============================================================================
//this function "envolope_angle" creates the maximal envoloppe of the angle 
//==============================================================================
function [Angle_envolope]=envolope_angle (Angle_measured)
    j=2;
    Angle_envolope = zeros(1,length(Angle_measured));
    Angle_envolope(1) =   Angle_measured(1);                          
    for i=2:length(Angle_measured)
        if(Angle_measured(i) > Angle_envolope(j-1)) then
            Angle_envolope(j) = Angle_measured(i);
            j=j+1;
        else
            Angle_envolope(j) =Angle_envolope(j-1);
            j=j+1;
        end
    end
endfunction
//==============================================================================
//this function "envolope_torque" creates the maximal envoloppe of the torque angle  
//==============================================================================
function [Torque_envolope]=envolope_torque(Torque_measured)
    j=2;
    Torque_envolope = zeros(1,length(Torque_measured)); 
    Torque_envolope(1) = Torque_measured(1);
    for i=2:length(Torque_measured)
        if( Torque_measured(i) > Torque_envolope(j-1)) then
            Torque_envolope(j) = Torque_measured(i);
            j=j+1;
        else
            Torque_envolope(j) = Torque_envolope(j-1);
            j=j+1;
        end
    end
endfunction
//==============================================================================
//this function "filtrage"" filters the data angle and torque 
//==============================================================================
function [Angle_filter,Torque_filter]=filtrage(Angle_envolope,Torque_envolope,Coefficient_filtrage)
    for i=1:length(Angle_envolope)
        if i==1 then 
            Torque_filter(i)=Torque_envolope(i);
            Angle_filter(i)=Angle_envolope(i);
        else

            Torque_filter(i)=Torque_filter(i-1)+Coefficient_filtrage*(Torque_envolope(i)-Torque_filter(i-1));
            Angle_filter(i)=Angle_filter(i-1)+Coefficient_filtrage*(Angle_envolope(i)-Angle_filter(i-1)); 
        end       
    end
endfunction
//==============================================================================
//the main task of this function "points_domaine_elastique" is select all points 
// of the electic domaine
//==============================================================================
function [Angle,Torque] =points_domaine_elastique(Angle_measured,Torque_measured,Pulse_torque,Target_torque)
    j=1; 
    for i=2: length(Angle_measured)
        if Pulse_torque(i)<> Pulse_torque(i-1) then 
            Angle(j)=Angle_measured(i-1);
            Torque(j)=Torque_measured(i-1);
            if Torque(j)>=Target_torque then
                break;
            end 
            j=j+1;
        end
    end
endfunction
//==============================================================================
//the main task of this function "points_previssage" is select all points 
// of the rundown domaine
//==============================================================================
function [Angle,Torque]=points_previssage (Angle_measured,Torque_measured,Pulse_torque)
    j=2;
    compteur=0;
    Angle(1)=Angle_measured(1);
    Torque(1)=Torque_measured(1);
    for i=2:length(Angle_measured)
        Angle(j)=Angle_measured(i);
        Torque(j)=Torque_measured(i);
        if Pulse_torque(i) <> Pulse_torque(i-1) then 
            compteur=compteur+1;
        end
        if compteur==1 then 
            break;
        end

        j=j+1;
    end
endfunction

//==============================================================================
//the main task of this function "pt_domain_elastic_tool_1" is select all points 
// of the elastic domaine using the points getted in the CV MONITOR 
//==============================================================================
function [Angle,Torque]=pt_domain_elastic_tool_1(Angle_measured,Torque_measured)
    j=1; 
    for i=1:length(Torque_measured)
        if Torque_measured(i)<>0 & Angle_measured(i)<>0 then 
            Angle(j)=Angle_measured(i);
            Torque(j)= Torque_measured(i);
            j=j+1;
        end

    end
endfunction
//==============================================================================
//the main task of this function "envolope_data_comtool" creates the envoloppe
//of the data getted from commtools
//==============================================================================
function [Angle,Torque,countt,verf]=envolope_data_comtool(Angle_measured,Pulse_torque,Torque_measured,Speed_measured)
    compteur=0;
    j=1;
countt=0;
verf=1;
for i=1:length(Pulse_torque)

    if i>1  then
        if Pulse_torque(i)<>Pulse_torque(i-1) then 
            Angle(j)=Angle_measured(i);
            Torque(j)=Pulse_torque(i);
            j=j+1
            compteur=compteur+1;
        end
        if compteur==0 then 
            if Speed_measured(i)<0 & Speed_measured(i-1)>0 then 
                Angle(j)=Angle_measured(i-1);
                Torque(j)=Torque_measured(i-1);
                j=j+1;
                countt=countt+1;
            end
            if countt==0 then 
                Angle(j)=Angle_measured(i-1);
                Torque(j)=Torque_measured(i-1);
                j=j+1 
            else
                Angle(j)=Angle(j-1);
                Torque(j)=Torque(j-1);
            end
        else
            Angle(j)=Angle(j-1);
            Torque(j)=Pulse_torque(i)
            j=j+1;
        end
    end
end

 

endfunction
//==============================================================================
//the main task of this function "find_position" is to place each variable 
//in the right position
//==============================================================================
function [position1,position2,position3,position4]= find_position(v)
    count=0
    for i=1:length(v) 
        if v(i)==10 then 
            position1=i+1;
            count=count+1;
        else
            if v(i)==14 then 
                position2=i+1;
                count=count+1; 
            else
                if v(i)==29 then
                    position3=i+1;
                    count=count+1;
                else
                    if v(i)==54 then 
                        position4=i+1;
                        count=count+1;
                    end


                end
            end
        end
    end
    if count==0  then
        messagebox('Data must countain Speed , Mecanical Angle  pulse torque and torque  ')
    else
        if count<4 then 
            messagebox("check data source Torque or Angle or Speed or pulse torque is missed ")
        end
    end
endfunction

//==============================================================================
//the main task of this function "find_puls_theshold_COM" is to find 
// the torque pulse threshold of data commtool
//==============================================================================
function [Tq_Puls_threshold,Ag_Puls_threshold_1, Ag_Puls_threshold_2]=find_puls_theshold_COM(Angle_measured,Pulse_torque)
    j=1;
    compteur=0;
    for i=2:length(Pulse_torque)

        if Pulse_torque(i)<>Pulse_torque(i-1) then 
            compteur=compteur+1;
        end
        if compteur==1 then 
            break;
        end

        j=j+1;
    end

    Tq_Puls_threshold = Pulse_torque(i);
    Ag_Puls_threshold_1=0;
    Ag_Puls_threshold_2=Angle_measured(i);

endfunction
//==============================================================================
//the main task of this function "find_puls_theshold_CVM" is to find 
// the torque pulse threshold of data CV MONITOR
//==============================================================================
function [Tq_Puls_threshold,Ag_Puls_threshold]=find_puls_theshold_CVM(Angle_measured,Torque_measured)
    compteur=0;
    for i=1:length(Angle_measured)

        if Angle_measured(i)<> 0 then 
            compteur=compteur+1;
        end
        if compteur==1 then 
            Tq_Puls_threshold=Torque_measured(i);
            Ag_Puls_threshold= Angle_measured(i);
            break;
        end
    end
endfunction
//==============================================================================
//the main task of this function "find_target_torque_COM" is to find 
// the TARGET TORQUE of data commtool
//==============================================================================
function [Target_torque,Target_Angle]=find_target_torque_COM(Angle_measured,Torque_measured)
    Target_torque=Torque_measured(1);

    for i=1:length(Torque_measured)

        if Torque_measured(i)>=Target_torque then 
            Target_torque=Torque_measured(i);
            Target_Angle=Angle_measured(i);
        end 



    end
endfunction
//==============================================================================
//the main task of this function "find_target_torque_CVM" is to find 
// the TARGET TORQUE of data CV monitor
//==============================================================================
function [Target_torque,Target_Angle]=find_target_torque_CVM(Angle_measured,Torque_measured)
    Target_torque=0;

    for i=1:length(Angle_measured)
        if Torque_measured(i)>Target_torque then 
            Target_torque=Torque_measured(i);
            Target_Angle=Angle_measured(i);
        end 
    end
endfunction


//==============================================================================
//the main task of this function "Angle_sign" is to change the sign of the angle  
// if it is negative 
//==============================================================================
function [Angle]=Angle_sign(Angle_measured)
    j=1; 
    indice_3=0;

    for i=1:length(Angle_measured)
        if Angle_measured(i)>0  then 
            indice_1=i; 
        else
            if Angle_measured(i)<-1  then
                indice_2=i;
                indice_3=indice_3+1;
            end
        end
        if indice_3==1 then 
            Delta=Angle_measured(indice_1)- Angle_measured(indice_2);
            break;
        end
    end 

    for i=1: length(Angle_measured)
        if indice_3==0 then 
            Angle(j)=Angle_measured(i)-Angle_measured(1);
            j=j+1
        else 
            if indice_3==1 then 
                if i<= indice_1 then 
                    Angle(j)=Angle_measured(i)-Angle_measured(1);
                    j=j+1;
                else
                    if i>=indice_2 then 
                        Angle(j)=Angle_measured(i)+Delta-Angle_measured(1);
                        j=j+1;
                    end
                end
            end
        end

    end
    Angle=Angle-Angle(1);
endfunction
//==============================================================================
//the main task of this function "sign_speed_randown" is change the speed sign 
// if it si negative at the beginning of the rundown domaine
//==============================================================================
function [Speed]=sign_speed_randown(Speed_measured)
    j=1;
    for i=1:length(Speed_measured)
        if i<=100 then 
            if Speed_measured(i)<0 then 
                Speed(j)=-Speed_measured(i);
                j=j+1;
            else
                Speed(j)=Speed_measured(i);
                j=j+1;
            end
        else
            Speed(j)=Speed_measured(i);
            j=j+1;
        end

    end
endfunction


//==============================================================================
//the main task of this function "select_point" is to select the points between
//50% and 100ù of target toque in ordre to calculate assembly stiffeness. 
//==============================================================================
function [Angle,Torque]=select_point(Anglein,Torquein,Target_torque)
    j=1;
    for i=1 : length(Torquein)
        if Torquein(i)>=0.5*Target_torque & Torquein(i)<=0.9*Target_torque then 
            Angle(j)=Anglein(i);
            Torque(j)=Torquein(i);
            j=j+1;
        end


    end
endfunction
//==============================================================================
//the main task of this function "assembly_type" is determine the type of the 
//Assembly
//==============================================================================
function [assembly_1,assembly_2]=assembly_type(x_elastic_origine,Target_angle,angle_median)

    if Target_angle > 360   then
        ecart=Target_angle -360;
        Target_angle=360;
        if x_elastic_origine> ecart
            x_elastic_origine=x_elastic_origine-ecart;
        end
        if angle_median> ecart then
            angle_median=angle_median-ecart;
        end
        difference_1= Target_angle-x_elastic_origine;
        difference_2= Target_angle-angle_median;
    else
        difference_1= Target_angle-x_elastic_origine;
        difference_2= Target_angle-angle_median;
    end

    if  difference_1 <= 30 then
        TmpMsg = msprintf('the assembly type using the ISO5393 norm is HARD ');
        Message_console(handles,TmpMsg);
        assembly_1=1;
    else
        if difference_1>30 & difference_1<=360 then 
            TmpMsg = msprintf('the assembly type using the ISO5393 norm is SOFT ');
            Message_console(handles,TmpMsg);
            assembly_1=2;
        end
    end
    if difference_2 <=30  then
        TmpMsg = msprintf('the assembly type using the VDI/VDE norm is HARD '); 
        Message_console(handles,TmpMsg);
        assembly_2=1;
    else
        if difference_2 >30 &  difference_2<=360 then 
            TmpMsg = msprintf('the assembly type using the VDI/VDE norm is SOFT ');
            Message_console(handles,TmpMsg);
            assembly_2=2;
        end

    end
endfunction
//==============================================================================
//the main task of this function "14" is to select the points 
// of the discontinuous domaine
//==============================================================================
function [Angle,Torque,Speed]=measure_domaine_puls(Angle_measured,Torque_measured,Pulse_torque,Speed_measured)
    count=0;
    j=1;
    for i=2:length(Pulse_torque)
        if Pulse_torque(i) <> Pulse_torque(i-1) then 
            count=count+1;
        end 
        if count <> 0 then 
            Angle(j)=Angle_measured(i);
            Torque(j)=Torque_measured(i);
            Speed(j)=Speed_measured(i);
            j=j+1;
        end


    end
    Angle=Angle-Angle(1); 

endfunction
//==============================================================================
//the main task of this function "simplified_representation" is to represent  
// simply the assembly with its remarkable points.
//==============================================================================
function [Angleout,Torqueout]=simplified_curve(Angle_measured,Torque_measured,Speed)
    compteur=0;
    k=2;

    Angleout(1)=Angle_measured(1);
    Torqueout(1)=Torque_measured(1);
    for i=1:length(Angle_measured)
        if i>1 then 
            if Speed(i)<0 & Speed(i-1)>0 then 
                compteur=compteur+1;
            else
                if compteur ==0 then 
                    Angleout(k)=Angle_measured(i);
                    Torqueout(k)=Torque_measured(i);
                    k=k+1;
                end
            end
        end               
    end

endfunction
//==============================================================================
//the main task of this function "filtrage_torque_rate" is to filter the 
// calculated torque rate. 
//==============================================================================
function [Torque_rate_out]=filtrage_torque_rate (Torque_rate,Coefficient_filtrage)
    for i=1:length(Torque_rate)
        if i==1 then
            Torque_rate_out(i)=Torque_rate(i);
        else 
            Torque_rate_out(i)= Torque_rate_out(i-1)+Coefficient_filtrage*(Torque_rate(i)- Torque_rate_out(i-1));
        end       
    end
endfunction
//==============================================================================
//the main task of this function "filtrage_torque_rate" is to filter the 
// calculated torque rate. 
//==============================================================================
function [TR_carre_cal,  coe_corr_lineaire]= Torquerate_moindre_carre(Angle,Torque,Angle_span,Taille_du_Buffer)

    indice_last_angle_step=1;// ns permet de conserver l'indice du i où moment on a arrivé à l'angle step
    Increment_buffer=1;
    Indice_first_buffer=1;
    Increment_Torque_rate=1;
    for i=1:length(Angle ) 
        // s'il s'agit du premiere indice on le stoke dans notre buffer 
        if i==1 then 
            Buffer_angle(Increment_buffer)=Angle (i)-Angle (1);
            Buffer_torque(Increment_buffer)=Torque (i);
            Indice_end_buffer =Increment_buffer;
            Increment_buffer=Increment_buffer+1;
            // sinon on paase au suivant 
        else 
            // on teste que la premiere difference entre le premier élement et celui de l'étape i est supérieur ou égale à l'angle step si ou onstoke l'element dans le buffer et on garde l'indice de l'élemnt dans la memoire

            if Angle(i)- Angle(indice_last_angle_step)>= Angle_span /(Taille_du_Buffer) // est que la difference est superieur ou égale à l'angle step
                Buffer_angle(Increment_buffer)=Angle (i)-Angle (1);                  //stokage de'angle dans le buffer 
                Buffer_torque(Increment_buffer)=Torque (i); 
                // tabtorqurate(Increment_buffer)  =Torque_rate_measured(i);                                //stokage du couple correspondant 
                indice_last_angle_step =i;                                                                // stokage de l'indice de l'élement qui correspond à l'angle step
                Indice_end_buffer=Increment_buffer;                                                               // stokage d'indice de dernier élement stoké dans le tableau
                Increment_buffer=Increment_buffer+1;
            end 

            if Buffer_angle(Indice_end_buffer)-Buffer_angle(Indice_first_buffer)>= Angle_span then                         // est ce que la difference est superieur ou égale à l'angle span
                // calcule du torque rate 
                cov_A_T =mean(Buffer_angle(Indice_first_buffer:Indice_end_buffer).*Buffer_torque(Indice_first_buffer:Indice_end_buffer))- mean(Buffer_angle(Indice_first_buffer:Indice_end_buffer))*mean(Buffer_torque(Indice_first_buffer:Indice_end_buffer));                                                             
                A=mean(Buffer_angle(Indice_first_buffer:Indice_end_buffer).^2)-(mean(Buffer_angle(Indice_first_buffer:Indice_end_buffer))).^2;
                T=mean(Buffer_torque(Indice_first_buffer:Indice_end_buffer).^2)-(mean(Buffer_torque(Indice_first_buffer:Indice_end_buffer))).^2;

                if A <> 0 then 
                    TR_carre_cal(Increment_Torque_rate)=  cov_A_T/ A;
                    coe_corr_lineaire(Increment_Torque_rate)=cov_A_T/sqrt(A*T);

                    Increment_Torque_rate=Increment_Torque_rate+1;                                                                 // indice d'incrementation du tableau torque rate 
                else 
                    TR_carre_cal(Increment_Torque_rate)= 0;
                    Increment_Torque_rate=Increment_Torque_rate+1; 
                end 
                Indice_first_buffer = Indice_first_buffer+1;             // déplacement au deuxieme élement du buffer      

            end 


        end
    end 
endfunction 
//==============================================================================
// this function "number_of_pulse" determines the number of pulses of tightening 
//==============================================================================
function [pulse_number]=number_of_pulse(Pulse_torque)
    pulse_number=0;
    for i=2:length(Pulse_torque)
        if Pulse_torque(i)<>Pulse_torque(i-1) then
            pulse_number=pulse_number+1;
        end
    end
    TmpMsg = msprintf('The number of pulse of tightening is =[%f]',pulse_number);
    Message_console(handles,TmpMsg);
endfunction
//==============================================================================
//this function "tightening_recommend" gives the user some advices that help him to improve 
// the capability of the tightening using the pulse tools. 
//=============================================================================
function [advice]=tightening_recommend(pulse_number, Pulse_Amplitude)
    if pulse_number==10 then
        TmpMsg = msprintf('tightening the assembly using this pulse tool respect the norm');
        Message_console(handles,TmpMsg);
        advice=1
    else
        if pulse_number > 10 then 
            if pulse_number <= 15 then 
                TmpMsg = msprintf('tightening the assembly using this pulse tool respect the norm');
                Message_console(handles,TmpMsg);
                advice=2
            else
                if Pulse_Amplitude == 100 then 
                    TmpMsg = msprintf('try to tight the assembly using another pulse tool ');
                    Message_console(handles,TmpMsg);
                    advice=3
                end
            end
        else
            if Pulse_Amplitude <= 100 & Pulse_Amplitude > 10
                TmpMsg = msprintf('decrease the amplitude pulse by 10 percent ');
                Message_console(handles,TmpMsg);
                advice=4

            end
        end

    end
endfunction
//==============================================================================
// this function ""calculates the error between the theoretical stiffness and 
// the calculated stiffeness
//==============================================================================
function [error_stiffness]=error_of_stifness(Raideur)
    [Target_torque,Target_Angle]=find_target_torque_COM(Angle_measured,Torque_measured)
    [Tq_Puls_threshold,Ag_Puls_threshold_1, Ag_Puls_threshold_2]=find_puls_theshold_COM(Angle_measured,Pulse_torque)
    Raideur_theo=(Target_torque-Tq_Puls_threshold)/(Target_Angle-Ag_Puls_threshold_2);
    if Raideur_theo > Raideur then
        error_stiffness= 100*(Raideur_theo-Raideur)/Raideur_theo;
        TmpMsg = msprintf('the error of calculation of the assembly stiffness is %f  ',error_stiffness);
        Message_console(handles,TmpMsg);
    else 
        error_stiffness= 100*(Raideur-Raideur_theo)/Raideur;
        TmpMsg = msprintf('the error of calculation of the assembly stiffness is %f ',error_stiffness);
        Message_console(handles,TmpMsg);
    end
endfunction

//==============================================================================
// the main task of this function "assembly_two_stiffeness" is to select the points
// by taking in account the variation of the stiffeness of the assembly
//==============================================================================
function [Angle,Torque]=pulse_domaine(Angle_measured,Pulse_torque)
    compteur=0;
    j=1;
    for i=2:length(Pulse_torque)
        if Pulse_torque(i) <> Pulse_torque(i-1) then 
            Angle(j)=Angle_measured(i);
            Torque(j)=Pulse_torque(i);
            compteur=compteur+1; 
            j=j+1;
        end
        if compteur<>0 then 
            Angle(j)=Angle(j-1);
            Torque(j)=Pulse_torque(i);
            j=j+1;
        end

    end
    Angle=Angle-Angle(1);
endfunction




//==============================================================================
// Définition de l'IHM
//=============================================================================
compteur=0;
while compteur==0
    TableMessage = ["  >","  >","  >","  >","  >","  >","  >","  >","  >"];
    f=figure('figure_position',[30,30],'figure_size',[1936,1096],'auto_resize','on','background',[33],'figure_name','Figure n°%d');
    //
    delmenu(f.figure_id,gettext('File'));
    delmenu(f.figure_id,gettext('?'));
    delmenu(f.figure_id,gettext('Tools'));
    //toolbar(f.figure_id,'off')
    //-----------------------------------------------------------------------------

    //-----------------------------------------------------------------------------
    handles.dummy = 0;
    handles.graph_angulaire= newaxes();
    handles.graph_angulaire.margins = [ 0 0 0 0];
    handles.graph_angulaire.axes_bounds = [0.09,0.04,0.35,0.35];
    title("Angular graph estimated")
    xlabel("Angle (°)")
    ylabel("Torque (N.m)")
    //-----------------------------------------------------------------------------
    handles.graph_torque_rate= newaxes();
    handles.graph_torque_rate.margins = [ 0 0 0 0];
    handles.graph_torque_rate.axes_bounds = [0.50,0.04,0.35,0.35];
    title("Torque rate graph");
    ylabel("Torque rate")
    //----------------------------------------------------------------------------
    handles.graph_simple_pres= newaxes();
    handles.graph_simple_pres.margins = [ 0 0 0 0];
    handles.graph_simple_pres.axes_bounds = [0.50,0.46,0.35,0.35];
    title("Simplified representation");
    xlabel(" Angle (°)")
    ylabel("Torque (N.m)")
    //-----------------------------------------------------------------------------
    handles.graph_tightening= newaxes();
    handles.graph_tightening.margins = [ 0 0 0 0];
    handles.graph_tightening.axes_bounds = [0.09,0.46,0.35,0.35];
    title("comparaison between measured and estimated curve")
    xlabel("Angle (°)")
    ylabel("Torque (N.m)")

    //-----------------------------------------------------------------------------

    // creation des bouton 
    handles.test=uicontrol(f,'unit','normalized','BackgroundColor',[-1,-1,-1],'Enable','on','FontAngle','normal','FontName','Tahoma','FontSize',[12],'FontUnits','points','FontWeight','normal','ForegroundColor',[-1,-1,-1],'HorizontalAlignment','center','ListboxTop',[],'Max',[1],'Min',[0],'Position',[0.09,0.11,0.1,0.05],'Relief','default','SliderStep',[0.01,0.1],'String','Run ','Style','pushbutton','Value',[0],'VerticalAlignment','middle','Visible','on','Tag','test','Callback','test_callback(handles)');
    handles.abort=uicontrol(f,'unit','normalized','BackgroundColor',[-1,-1,-1],'Enable','on','FontAngle','normal','FontName','Tahoma','FontSize',[12],'FontUnits','points','FontWeight','normal','ForegroundColor',[-1,-1,-1],'HorizontalAlignment','center','ListboxTop',[],'Max',[1],'Min',[0],'Position',[0.20,0.11,0.1,0.05],'Relief','default','SliderStep',[0.01,0.1],'String','Quit','Style','pushbutton','Value',[0],'VerticalAlignment','middle','Visible','on','Tag','test','Callback','abort_callback(handles)');
    handles.actualiser=uicontrol(f,'unit','normalized','BackgroundColor',[-1,-1,-1],'Enable','on','FontAngle','normal','FontName','Tahoma','FontSize',[12],'FontUnits','points','FontWeight','normal','ForegroundColor',[-1,-1,-1],'HorizontalAlignment','center','ListboxTop',[],'Max',[1],'Min',[0],'Position',[0.31,0.11,0.1,0.05],'Relief','default','SliderStep',[0.01,0.1],'String','Update','Style','pushbutton','Value',[0],'VerticalAlignment','middle','Visible','on','Tag','test','Callback','Actualisation_callback');
    //------------------------------------------------------------------------------
    handles.Run_down=uicontrol(f,"style","radiobutton","string","Torque rate rundow domaine",'Position',[115 45 200 15],'groupname','grouptest');
    handles.Discountinuous_tightenin=uicontrol(f,"style","radiobutton","string","Torque rate discontinous domain",'Position',[365 45 200 15],'groupname','grouptest');
    handles.elastic_domaine=uicontrol(f,"style","radiobutton","string","Torque rate elastique domain",'position',[615 45 200 15],'groupname','grouptest')
    //------------------------------------------------------------------------------


    //------------------------------------------------------------------------------
    graph1=handles.graph_angulaire;
    graph2=handles.graph_torque_rate;
    graph3=handles.graph_tightening;
    graph4=handles.graph_simple_pres

    //------------------------------------------------------------------------------
    check_Run_down=handles.Run_down;
    check_Disc_tightening=handles.Discountinuous_tightenin;
    check_elastic_domaine=handles.elastic_domaine;

    //------------------------------------------------------------------------------

    set(check_Run_down,"Enable","on","value","0"); 
    set(check_Disc_tightening,"Enable","on","value","0");
    set(check_elastic_domaine,"Enable","on","value","0");

    //------------------------------------------------------------------------------
    w = createWindow()
    w.axes_size = [800 800];
    handles.message_consol=uicontrol(w,"style","text","string","",'Position',[10 380 800 20],'fontsize',14,'FontWeight','demi');

    TmpMsg = msprintf('Select analysis options \n');
    Message_console(handles,TmpMsg);
    TmpMsg = msprintf('Start analysis to choose a file \n');
    Message_console(handles,TmpMsg);
    //============================================================================//
    //                 PROGRAMME PRINCIPALE                                       //
    //===========================================================================//    

    Run_Analyze = 0;
    while Run_Analyze == 0
        sleep(100);


    end
    //==============================================================================

    if (Run_Analyze == 1) then

        a=sca(graph1);
        delete(a.children);
        a=sca(graph2);
        delete(a.children);
        a=sca(graph1);

        //==============================================================================
        // chargement du fichier 
        //==============================================================================
        text_file=uigetfile(["*.csv";"*.txt"]);
        fd=mopen(text_file); //ouverture du fichier

        he=mgetl(fd,1); // récuperation de la 1ere ligne
        mclose(fd); // fermeture du fichier
        l1=he(1);
        words=strsplit(l1);
        words1=strcmp(words,'C');
        //==============================================================================
        //                          type du file
        //==============================================================================        
        if l1=='sep=;' then
            file_type=1;
        else 
            if words1(1)==0 then 
                file_type=2;  

            else 
                file_type=0;

            end
        end          
        //------------------------------------------------------------------------------
        TmpMsg = msprintf('file type : %s',string( file_type));
        Message_console(handles,TmpMsg);
        //==============================================================================
        //             récupération des données de Comm Tool
        //==============================================================================
        if  file_type==2 then // récupération des données de Comm Tool
            TmpMsg = msprintf('getting data from Comm Tool');
            Message_console(handles,TmpMsg);
            labels=["Angle span (degre)","length of Buffer","Coefficient of filter (<1)","Torque multiplication factor","Ratio","Speed multiplication factor","Pulse Amplitude (<=100%) "];
            [ok,Angle_span,Taille_du_Buffer,Coefficient_filtrage,Torque_factor,Ratio,Speed_factor, Pulse_Amplitude ]= getvalue("Parameters",labels,list("vec",1,"vec",1,"vec",1,"vec",1,"vec",1,"vec",1,"vec",1),["5","99","0.1","15.02/821","8.5","0.115","100"])
            fd = mopen(text_file); //ouverture du fichier
            he = mgetl(fd,9); // récuperation du header
            l7 = he(7); // récupération de la ligne 7
            TimeBase = strtod(part(l7,strindex(l7,':')+1:strindex(l7,':')+5)); // récupération de la fréquence d'échantillonnage
            l6 = he(6); // récupération de la ligne 6
            Nb_points=strtod(part(l6,strindex(l6,':')+1:strindex(l7,':')+8)); // récupération du nombre de points
            frequence = 1000000/TimeBase;
            mclose(fd); // fermeture du fichier
            header = 9; // 9 lignes de commentaires avant les données dans le fichier de exporté de Comm Tool
            TablPoints=csvRead(text_file, ";",'double',[],[],'/Data RAW/',[],header); //récupération du fichier de points de Comm Tool
            //==============================================================================
            ligne_1=tokens(he(1),':');
            ligne_2=tokens(he(2),':');
            ligne_3=tokens(he(3),':');
            ligne_4=tokens(he(4),':');
            v=[ligne_1(2),ligne_2(2),ligne_3(2),ligne_4(2)]
            v=strtod(v);
            [position1,position2,position3,position4]=find_position(v)
            time     = TablPoints(1:Nb_points,1)/1000; // définition du vecteur temps
            length_t = size(time,'*');
            Torque_measured  = Torque_factor*TablPoints(1:Nb_points,position2); // vecteur du couple Torque
            Angle_measured1    =TablPoints(1:Nb_points,position3)*360/(Ratio*4096); // vecteur de l'angle Mecanical turns
            [Angle_measured]=Angle_sign(Angle_measured1);
            Speed_measured= Speed_factor*TablPoints(1:Nb_points,position1);// vecteur de vitesse 
            Pulse_torque=Torque_factor*TablPoints(1:Nb_points,position4);
            //==============================================================================
            //                 récuperation des données CVI MONITOR
            //==============================================================================
        elseif  file_type == 1 then 
            TmpMsg = msprintf('getting data from CVI MONITOR');
            Message_console(handles,TmpMsg);

            fd = mopen(text_file); // open the file 

            he = mgetl(fd,3);

            mclose(fd); // close the file
            header = 3; //  get all lines before data in the file CV monitor 
            TablPoints=csvRead(text_file, ";",'string','string',[],'/Data RAW/',[],header); // get the data from file CV MONITOR 
            strsubst( text_file, ",", ".");
            TablPoints_size=size(TablPoints);


            time     = TablPoints(:,1) // definethe time vector 
            length_t = size(time,'*');
            Angle_measured1    = TablPoints(1:TablPoints_size(1),2) ; // define the angle vector
            Angle_measured=strtod( Angle_measured1)
            Torque_measured1   = TablPoints(1:TablPoints_size(1),3); // define the torque vector 
            Torque_measured =strtod( Torque_measured1);
            Current_measured1 =TablPoints(1:TablPoints_size(1),4);// define the electric currrent vector 
            Current_measured =strtod(  Current_measured1);
            Speed_measured1= TablPoints(1:TablPoints_size(1),5);//define the rundown speed vector 
            Speed_measured =strtod(  Speed_measured1);
        end,


        //============================================================================//
        //          CALCULATION OF THE TORQUE RATE                                    //
        //============================================================================//

        if file_type == 2 then 

            //-----------------------------------------------------------------------------
            if check_Run_down.value == 1  then
                //--------------------------------------------------------------------------
                [Angle,Torque]=points_previssage (Angle_measured,Torque_measured,Pulse_torque); // select hte points of run-down
                //--------------------------------------------------------------------------
                [Torque_rate]=Torquerate_method_median(Angle,Torque,Angle_span,Taille_du_Buffer); // calculate  the torque rate of the ey
                //--------------------------------------------------------------------------
                a=sca(graph1);
                plot(Angle,Torque,"-.k");
                a=sca(graph2);
                plot(Torque_rate,"-.r");
                //--------------------------------------------------------------------------
                //clear Angle;
                //clear Torque;
                //------------------------------------------------------------------------------        
            end

            if  check_elastic_domaine.value==1 then 

                //--------------------------------------------------------------------------
                [Target_torque,Target_Angle]=find_target_torque_COM(Angle_measured,Torque_measured)
                //--------------------------------------------------------------------------
                [Angle_elastic,Torque_elastic] =points_domaine_elastique(Angle_measured,Torque_measured,Pulse_torque,Target_torque); 
                //--------------------------------------------------------------------------
                [Torque_rate]=Torque_rate_discontinu(Angle_elastic,Torque_elastic);
                //--------------------------------------------------------------------------
                [ count, median_angle_1,median_torque_1, median_angle_2,median_torque_2]=poly_stiffness(Angle_elastic,Torque_elastic,Torque_rate);
                [Angle1,Torque1]=select_point(Angle_elastic,Torque_elastic,Target_torque)
                //--------------------------------------------------------------------------
                if count==0 then
                [Raideur]=Torque_rate_discontinu(Angle1,Torque1);
                ///calculation of the stiffeness of assembly =======================================
                raideur_Assemblage=abs(mean(Raideur));
                TmpMsg = msprintf('the stiffness of the Assembly is :  --> .%f ',raideur_Assemblage);
                Message_console(handles,TmpMsg)
                end
                //==============================================================
                a=sca(graph2);
                plot(Torque_rate,"-.k");
                a=sca(graph1);
                plot(Angle_elastic,Torque_elastic,"-.r");
                clear Angle;
                clear Torque;
            end

            if  check_Disc_tightening.value==1 then 
                [Speed]=sign_speed_randown(Speed_measured) // get red of the negative speed during the ryn 
                //--------------------------------------------------------------------------
                [Angle_real,Torque_real,Speed_real]=measure_domaine_puls(Angle_measured,Torque_measured,Pulse_torque,Speed) // get discontinuous tightening data 
                //--------------------------------------------------------------------------
                [Target_torque,Target_Angle]=find_target_torque_COM(Angle_real,Torque_real) // find target torqur and targrt angle 
                //--------------------------------------------------------------------------
                [Anglefr, Torquefr, countt,verf]=envolope_data_comtool(Angle_measured,Pulse_torque,Torque_measured,Speed) // get the envollope of the data 
                //--------------------------------------------------------------------------
                [Angle_filter,Torque_filter]=filtrage(Anglefr,Torquefr,Coefficient_filtrage);//
                //--------------------------------------------------------------------------
                [Torque_rate]=Torquerate_method_median(Angle_filter,Torque_filter,Angle_span,Taille_du_Buffer)
                //-------------------------------------------------------------------------------------
                 [Torque_rate1]=Torquerate_method_median(Anglefr,Torquefr,Angle_span,Taille_du_Buffer)
                 //----------------------------------------------------------------------------------
                [TR_carre_cal,  coe_corr_lineaire]= Torquerate_moindre_carre(Angle_filter,Torque_filter,Angle_span,Taille_du_Buffer)
                //----------------------------------------------------------------------------------
                [Angle_elastic,Torque_elastic] =points_domaine_elastique(Angle_measured,Torque_measured,Pulse_torque,Target_torque); 
                //------------------------------------------------------------------
                [Angle1,Torque1]=select_point(Angle_elastic,Torque_elastic,Target_torque)
                //-------------------------------------------------------------------
                [Torque_rate_1]=Torque_rate_discontinu(Angle_elastic,Torque_elastic);
                //-----------------------------------------------------------------------------------
                [ count, median_angle_1,median_torque_1, median_angle_2,median_torque_2]=poly_stiffness(Angle_elastic,Torque_elastic,Torque_rate_1);
                //--------------------------------------------------------------------------
                if count==0 then 
                    raideur_Assemblage=abs(mean(Torque_rate_1));
                    //-----------------------------------------------------------------------------------
                    TmpMsg = msprintf('the stiffness of the Assembly is :  --> .%f ',raideur_Assemblage);
                    Message_console(handles,TmpMsg);
                end
                //==========================================================
                a=sca(graph1);
                 plot(Anglefr,Torquefr,"k")
              //  plot(Angle_filter,Torque_filter,"-r")
                plot(Angle_measured,Torque_measured,"-g")
                legend("estimation whith speed sign","measured curve",[2])
                //==========================================================
                a=sca(graph2);
                plot(Torque_rate,"-k")
                plot(TR_carre_cal,"-r")
                plot(Torque_rate,"-b")
                legend("Torque rate median","TR less square",[2]);
                //==========================================================
                //clear Angle;
                //clear Torque;

            end
            //==================================================================

            [Speed]=sign_speed_randown(Speed_measured) // get red of the negative value of speed in run dow area 
            //------------------------------------------------------------------
            [Angle_real,Torque_real,Speed_real]=measure_domaine_puls(Angle_measured,Torque_measured,Pulse_torque,Speed)// select the data of discontinuous domaine  
            //------------------------------------------------------------------
            [Target_torque_1,Target_Angle_1]=find_target_torque_COM(Angle_real,Torque_real)// find target torque and target angle for discontinuous tightening data 
            //------------------------------------------------------------------
            [Target_torque_2,Target_Angle_2]=find_target_torque_COM(Angle_measured,Torque_measured)// find target torque and target angle 
            //------------------------------------------------------------------
            [Angleout,Torqueout]=simplified_curve(Angle_measured,Torque_measured,Speed)// get the data of run down area  
            //------------------------------------------------------------------
            [Tq_Puls_threshold,Ag_Puls_threshold_1, Ag_Puls_threshold_2]=find_puls_theshold_COM(Angle_measured,Pulse_torque) // find the puls torque 
            //--------------------------------------------------------------------------
            [Angle,Torque]=pulse_domaine(Angle_measured,Pulse_torque)
            //--------------------------------------------------------------------------
            [Angle_filter,Torque_filter]=filtrage(Angle,Torque,Coefficient_filtrage)
            //------------------------------------------------------------------
            [Angle_elastic,Torque_elastic] =points_domaine_elastique(Angle_measured,Torque_measured,Pulse_torque,Target_torque_1); 
            //------------------------------------------------------------------
            [Torque_rate_2]=Torque_rate_discontinu(Angle_elastic,Torque_elastic);
            //-----------------------------------------------------------------------------------
            [ count, median_angle_1,median_torque_1, median_angle_2,median_torque_2]=poly_stiffness(Angle_elastic,Torque_elastic,Torque_rate_2);
            //------------------------------------------------------------------
            if count ==0 then
                x_elastic_origin_2= Ag_Puls_threshold_2-Tq_Puls_threshold*((Ag_Puls_threshold_2-Target_Angle_2)/(Tq_Puls_threshold-Target_torque_2));
                if Ag_Puls_threshold_2 > x_elastic_origin_2 then 
                    x_elastic_line_2=[x_elastic_origin_2,Ag_Puls_threshold_2, Target_Angle_2]
                    y_elastic_line_2=[0,Tq_Puls_threshold,Target_torque_2]
                else 
                    x_elastic_line_2=[Ag_Puls_threshold_2,x_elastic_origin_2, Target_Angle_2]
                    y_elastic_line_2=[Tq_Puls_threshold,0,Target_torque_2]
                end
                Angle_present_2=[Angleout',Ag_Puls_threshold_2,Target_Angle_2]
                Torque_present_2=[Torqueout',Tq_Puls_threshold,Target_torque_2]
                //--------------------------------------------------------------------------
                x_elastic_line_1=[Ag_Puls_threshold_1, Target_Angle_1]
                y_elastic_line_1=[Tq_Puls_threshold,Target_torque_1]


                Angle_present_1=[Ag_Puls_threshold_1,Target_Angle_1]
                Torque_present_1=[Tq_Puls_threshold,Target_torque_1]

                //==================================================================
                a=sca(graph3)
                plot(Angle_filter,Torque_filter,"-.r")
                plot(x_elastic_line_1,y_elastic_line_1,"-.k")
                legend("estimated curve","theoritical curve",[2])
                a=gca();
                a.data_bounds=[0 0 ;Target_Angle_1+20  Target_torque_1+10 ];
                plot2d([Ag_Puls_threshold_1,Ag_Puls_threshold_1,1],[0,Tq_Puls_threshold,Tq_Puls_threshold],style=3)
                xstring(Ag_Puls_threshold_1+1,Tq_Puls_threshold,"(pulse threshold) ")
                TmpMsg = msprintf('The coordinates of the puls threshold point are : ( Angle puls threshold,Torque pulse threshold)=[%f° %f N.m ]',Ag_Puls_threshold_1,Tq_Puls_threshold);
                Message_console(handles,TmpMsg);
                plot2d([Target_Angle_1,Target_Angle_1,1],[0,Target_torque_1,Target_torque_1],style=3)
                xstring(Target_Angle_1,Target_torque_1,"(Target )")
                TmpMsg = msprintf('The coordinates of the target point are : ( Target Angle,Target Torque )=[%f° %f N.m ]',Target_Angle_1,Target_torque_1);
                Message_console(handles,TmpMsg);
                //==================================================================
                a=sca(graph4)
                plot(Angle_present_2,Torque_present_2,"-r")
                plot(x_elastic_line_2,y_elastic_line_2,"-k")
                //--------------------------------------------------------------------------
                a=gca();
                a.data_bounds=[0 0 ;Target_Angle_2+20  Target_torque_2+10 ];
                //--------------------------------------------------------------------------
                plot2d([Ag_Puls_threshold_2,Ag_Puls_threshold_2,1],[0,Tq_Puls_threshold,Tq_Puls_threshold],style=3)
                xstring(Ag_Puls_threshold_2-15,Tq_Puls_threshold,"( pulse threshold) ")
                //--------------------------------------------------------------------------
                plot2d([Target_Angle_2,Target_Angle_2,1],[0,Target_torque_2,Target_torque_2],style=3)
                xstring(Target_Angle_2-10,Target_torque_2,"(Target )")
                xstring(x_elastic_origin_2,0.25,"elastic origine")
                //--------------------------------------------------------------------------
                TmpMsg = msprintf('The coordinates of the elastic origin point are : ( elastic origin )=[%f° %f N.m ]',x_elastic_origin_2,0);
                Message_console(handles,TmpMsg);
                //------------------------ ASSEMBLY TYPE
                [Angle_out,Torque_out]=select_point(Angle,Torque,Target_torque_2)
                angle_median=Angle_out(1);
                [assembly_1,assembly_2]=assembly_type(x_elastic_origin_2,Target_Angle_2,angle_median)
                //--------------------------NUMBER OF PULSE OF TIGHTENING 
                [pulse_number]=number_of_pulse(Pulse_torque)
                //------------------------------ TIGHTENING RECOMMENDATION 
                [advice]=tightening_recommend(pulse_number, Pulse_Amplitude)
                if check_elastic_domaine.value==1 | check_Disc_tightening.value==1 then 
                    //------------------------------- error of stiffness
                    [error_stiffness]=error_of_stifness(raideur_Assemblage)
                end
            else
                if count==1 then 
                     //--------------------------NUMBER OF PULSE OF TIGHTENING 
                    [pulse_number]=number_of_pulse(Pulse_torque)
                    //------------------------------ TIGHTENING RECOMMENDATION 
                    [advice]=tightening_recommend(pulse_number, Pulse_Amplitude)
                    raideur_Assemblage1=(Tq_Puls_threshold-median_torque_1)/(Ag_Puls_threshold_2-median_angle_1)
                    //--------------------------------------------------------------------------
                    TmpMsg = msprintf('the stiffness of the first part of the Assembly is :  --> .%f ',raideur_Assemblage1);
                    Message_console(handles,TmpMsg);
                    raideur_Assemblage2= (Target_torque_2-median_torque_2)/(Target_Angle_2-median_angle_2);
                    //--------------------------------------------------------------------------
                    TmpMsg = msprintf('the stiffness of the second part of the Assembly is :  --> .%f ',raideur_Assemblage2);
                    Message_console(handles,TmpMsg);
                   
                   
                end
            end 

            //==============================================================================                                          
        else 
            //====================== calcule du torque rate pour les données CV MONITOR=====
            if file_type == 1 then 
                labels=["Angle span (degre)","length of Buffer","Coefficient of filter (<1)"];
                [ok,Angle_span,Taille_du_Buffer,Coefficient_filtrage]= getvalue("Parameters",labels,list("vec",1,"vec",1,"vec",1),["5","99","0.01"]);
                if check_Run_down.value == 1  then

                    messagebox('the angle equals to zero, impossible to compete the torque rate for this case ');
                end

                if check_elastic_domaine.value==1 then 
                    //--------------------------------------------------------------------------
                    [Angle,Torque]=pt_domain_elastic_tool_1(Angle_measured,Torque_measured)
                    //--------------------------------------------------------------------------
                    [Target_torque,Target_Angle]=find_target_torque_CVM(Angle_measured,Torque_measured) 
                    //--------------------------------------------------------------------------
                    [Torque_rate]=Torque_rate_discontinu(Angle,Torque); // appel de la fonction pour calculer le torque rate en previssage 
                    //--------------------------------------------------------------------------
                    [Angle2,Torque2]=select_point(Angle,Torque,Target_torque)// selection des points entre 50% et 95% du target torque 
                    //--------------------------------------------------------------------------
                    [Raideur]=Torque_rate_discontinu(Angle2,Torque2);//calcule de la raideur 
                    ///afficher le torque rate =====================================================
                    a=sca(graph2);
                    plot(Torque_rate,"-.r");
                    // affichage angulaire couple angle 
                    a=sca(graph1);
                    plot(Angle,Torque,"-.k");
                    ///calcule de la raideur de l'assemblage =======================================
                    raideur_Assemblage=abs(mean(Raideur));
                    TmpMsg = msprintf('the stiffness of the assembly is :  --> .%f ',raideur_Assemblage);
                    Message_console(handles,TmpMsg);
                    clear Angle;
                    clear Torque;
                end 
                if  check_Disc_tightening.value==1 then 

                    //-------------------------------------------------------------------------- 
                    [Angle_envolope]=envolope_angle (Angle_measured);
                    //--------------------------------------------------------------------------
                    [Torque_envolope]=envolope_torque(Torque_measured);
                    //--------------------------------------------------------------------------
                    [Target_torque,Target_Angle]=find_target_torque_CVM(Angle_measured,Torque_measured);
                    //--------------------------------------------------------------------------
                    [Angle_filter,Torque_filter]=filtrage(Angle_envolope,Torque_envolope,Coefficient_filtrage);
                    //--------------------------------------------------------------------------
                    [TR_carre_cal,  coe_corr_lineaire]= Torquerate_moindre_carre(Angle_filter,Torque_filter,Angle_span,Taille_du_Buffer)
                    //--------------------------------------------------------------------------
                    [Angle2,Torque2]=select_point(Angle_envolope,Torque_envolope,Target_torque)
                    //-------------------------------------------------------------------------- 
                    [Raideur]=Torque_rate_discontinu(Angle2,Torque2);//calcule de la raideur 
                    raideur_Assemblage=abs(mean(Raideur));
                    TmpMsg = msprintf('The stiffness of the Assembly is :  --> .%f ', raideur_Assemblage);
                    Message_console(handles,TmpMsg);
                    ///afficher le torque rate ====================================================
                    a=sca(graph2);
                    plot(TR_carre_cal,"-k.");
                    a=sca(graph1);
                    plot(Angle_filter,Torque_filter,"-r.");

                end 
                //--------------------------------------------------------------------------   
                if check_elastic_domaine.value==1 | check_Disc_tightening.value==1 then 
                    [Target_torque,Target_Angle]=find_target_torque_CVM(Angle_measured,Torque_measured) 
                    //--------------------------------------------------------------------------
                    [Tq_Puls_threshold,Ag_Puls_threshold]=find_puls_theshold_CVM(Angle_measured,Torque_measured)
                    //--------------------------------------------------------------------------
                    x_elastic_origin= Ag_Puls_threshold-Tq_Puls_threshold*((Ag_Puls_threshold-Target_Angle)/(Tq_Puls_threshold-Target_torque));
                    if Ag_Puls_threshold > x_elastic_origin then 
                        x_elastic_line=[x_elastic_origin,Ag_Puls_threshold, Target_Angle]
                        y_elastic_line=[0,Tq_Puls_threshold,Target_torque]
                    else 
                        x_elastic_line=[Ag_Puls_threshold,x_elastic_origin, Target_Angle]
                        y_elastic_line=[Tq_Puls_threshold,0,Target_torque]
                    end
                    //--------------------------------------------------------------------------
                    a=sca(graph3)

                    plot(x_elastic_line,y_elastic_line,"-r")
                    plot(Angle_filter,Torque_filter,"-k.");
                    a=gca();
                    a.data_bounds=[0 0 ;Target_Angle+10 Target_torque+10 ];
                    plot2d([Ag_Puls_threshold,Ag_Puls_threshold,1],[0,Tq_Puls_threshold,Tq_Puls_threshold],style=3)
                    xstring(0.1,Tq_Puls_threshold+0.1,"Torque pulse threshold")
                    xstring(Ag_Puls_threshold+0.1,0.25,"Angle pulse threshold")
                    TmpMsg = msprintf('The coordinates of the puls threshold point are : ( Angle puls threshold,Torque pulse threshold)=[%f° %f N.m ]',Ag_Puls_threshold,Tq_Puls_threshold);
                    Message_console(handles,TmpMsg);
                    plot2d([Target_Angle,Target_Angle,1],[0,Target_torque,Target_torque],style=3)
                    xstring(0.1,Target_torque+0.1,"Target Torque ")
                    xstring( Target_Angle+0.1,0.25,"Target Angle " )
                    TmpMsg = msprintf('The coordinates of the target point are : ( Target Angle,Target torque )=[%f° %f N.m ]',Target_Angle,Target_torque);
                    Message_console(handles,TmpMsg);
                    //--------------------------------------------------------------------------
                    a=sca(graph4)
                    plot(x_elastic_line,y_elastic_line,"-k")
                    a=gca();
                    a.data_bounds=[0 0 ;Target_Angle+20  Target_torque+10 ];
                    plot2d([Ag_Puls_threshold,Ag_Puls_threshold,1],[0,Tq_Puls_threshold,Tq_Puls_threshold],style=3)
                    xstring(Ag_Puls_threshold-15,Tq_Puls_threshold,"( pulse threshold) ")
                    plot2d([Target_Angle,Target_Angle,1],[0,Target_torque,Target_torque],style=3)
                    xstring(Target_Angle-10,Target_torque,"(Target )")
                    xstring(x_elastic_origin,0.25,"elastic origine")
                    TmpMsg = msprintf('The coordinates of the elastic origin point are : ( elastic origin )=[%f° %f N.m ]',x_elastic_origin,0);
                    Message_console(handles,TmpMsg);
                    //--------------------------------------------------------------------------
                    // assemby type 
                    [Angle_out,Torque_out]=select_point(Angle_filter,Torque_filter,Target_torque)
                    angle_median=Angle_out(1);
                    [assembly_1,assembly_2]=assembly_type(x_elastic_origin,Target_Angle,angle_median)
                end
            else 
                if file_type==0 then 
                    messagebox("unrecognized file, the data must be from COMMTOOL or CV-Monitor, press the bouton Quit and try again with other file");
                end 

            end
        end    
        //==============================================================================
        //=========== actualisation du programme =======================================
        Actualiser = 0;
        while Actualiser == 0
            sleep(100);
        end 
        //============================================================================== 
        if Actualiser==1 then 
            compteur=0;
        end
        close(f);
        close(w)  ;  



    end 
end

















