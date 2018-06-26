classdef Parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hCom
        delta
        omega
        singleSupportDuration
        doubleSupportDuration
        firstSupportFoot
        predictionTime
        nSamples
        footSize
        feasibilityL
        feasibilityX
        feasibilityY
        tailType
        
        boxA
        boxB
        boxC
        
        sampleScaling
        deltas
        
        %swing foot
        alpha
        beta
    end
    
    methods
        function obj = Parameters(hCom, delta, singleSupportDuration, doubleSupportDuration, firstSupportFoot, ...
                                  predictionTime, nSamples, footSize, feasibilityL, feasibilityX, feasibilityY, tailType)

            obj.hCom = hCom;
            obj.delta = delta;
        	obj.singleSupportDuration = singleSupportDuration;
            obj.doubleSupportDuration = doubleSupportDuration;
            obj.firstSupportFoot = firstSupportFoot;
            obj.predictionTime = predictionTime;
            obj.nSamples = nSamples;
            obj.footSize = footSize;
            obj.feasibilityL = feasibilityL;
            obj.feasibilityX = feasibilityX;
            obj.feasibilityY = feasibilityY;
            obj.tailType = tailType;
            
            obj.boxA = obj.feasibilityL - feasibilityX;
            obj.boxB = obj.feasibilityL + feasibilityX;
            obj.boxC = obj.feasibilityY;
            
            obj.omega = sqrt(9.81/obj.hCom);
            
            % Compute the scaling for the time discretization
            poly(1) = 1;
            poly(obj.nSamples) = -obj.predictionTime/obj.delta;
            poly(obj.nSamples+1) = obj.predictionTime/obj.delta - 1;
            q = roots(poly);
            q = q(imag(q)==0);
            q = max(q);
            obj.sampleScaling = q;
            obj.deltas = obj.delta*ones(obj.nSamples);%0.01*q.^(0:obj.nSamples-1);            
        end
    end
    
end

