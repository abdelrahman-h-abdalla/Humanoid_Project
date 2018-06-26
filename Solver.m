classdef Solver < dynamicprops
    
    properties(GetAccess=public)
        % gains
        QZd
        QJ
        QVx
        QVy
        QZx
        QZy
        Qfsx
        Qfsy
        QCPx
        QCPy
        
        % parameters
        ss_time
        ds_time
        S_samples
        D_samples
        M
        N_samples
        F_samples
        delta
        omega
        w
        preassigned_footsteps_matrix
        footsteps_are_preassigned
        first_support_foot
        box_a
        box_b
        box_c
        
        % state
        x
        xd
        y
        yd
        t_curr
        zx
        zy
        theta
        iter
        iterF
        fxc
        fyc
        footstep
        vref_x
        vref_y
        vref_omega
        
        %swing foot state
        x_l
        xd_l
        y_l
        yd_l
        zx_l
        zy_l
        
        x_r
        xd_r
        y_r
        yd_r
        zx_r
        zy_r
        
        % plot utilities
        zd
        pred_zx
        pred_zy
        pred_zx_body
        pred_zy_body
        pred_zx_foot
        pred_zy_foot
        pred_fs
        plot_limits
        plot_options
        
        % storage arrays
        x_store
        y_store
        xd_store
        yd_store
        zx_store
        zy_store
        zxd_store
        zyd_store
        footsteps
        
        % solving utilities
        HQ
        fQ
        Vu
        Vs
        Pu
        Ps
        Cc
        A_zmp
        b_zmp
        A_fs
        b_fs
        Aeq
        beq
        predicted_rotations
        
        % number of sampling instants
        W
        pre_time
        p
        P
        deltas
        time_discretization_mask
        q
        
        parameters
        
        current_support_foot
        L
        
        orientations
        plot_handles
    end
    
    methods
        function obj = Solver(state, parameters)

            obj.L = parameters.feasibilityL-parameters.feasibilityY;
            obj.box_a = obj.L-parameters.feasibilityY;
            obj.box_b = parameters.feasibilityL+parameters.feasibilityY;
            obj.box_c = parameters.feasibilityX;
            
            % Initialize global variables
            obj.delta = parameters.delta;
            obj.omega = parameters.omega;
            obj.w = parameters.footSize;
            obj.ss_time = parameters.singleSupportDuration;
            obj.ds_time = parameters.doubleSupportDuration;
            obj.pre_time = parameters.predictionTime;
            
            obj.parameters = parameters;
            
            % Initialize state
            state = num2cell(state');
            [obj.x, obj.xd, xdd, obj.y, obj.yd, ydd, obj.theta, obj.t_curr] = state{:};
            
            first_foot_position = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)]*[0;obj.L/2];
            
            obj.fxc = obj.x + first_foot_position(1);
            obj.fyc = obj.y + first_foot_position(2);
            
            obj.x = 0;
            obj.y = -obj.L/2;
            
            obj.zx = obj.x - (1/obj.omega^2)*xdd;
            obj.zy = obj.y - (1/obj.omega^2)*ydd;
            
            obj.N_samples = round(obj.pre_time/obj.delta);
            obj.S_samples = round(obj.ss_time/obj.delta);
            obj.D_samples = round(obj.ds_time/obj.delta);
            obj.F_samples = (obj.S_samples+obj.D_samples)*2;
            obj.M = ceil(obj.N_samples/(obj.S_samples+obj.D_samples));
            N=obj.N_samples;S=obj.S_samples;D=obj.D_samples;F=obj.F_samples;M=obj.M;
            
            % Variable discretization
            obj.W = parameters.nSamples;
            
            obj.deltas = parameters.deltas;
            
            obj.p = ones(obj.W,1);
            
            obj.P = [];
            for i = 1:obj.W
                obj.P(i,1:i) = obj.deltas(1:i);
            end
            
            % Initialize storage vectors           
            rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
            
            pos_abs = [obj.fxc;obj.fyc] + rot*[obj.x;obj.y];
            vel_abs = rot*[obj.xd;obj.yd];
            zmp_abs = [obj.fxc;obj.fyc] + rot*[obj.zx;obj.zy];
            obj.x_store(1) = pos_abs(1);
            obj.y_store(1) = pos_abs(2);
            obj.xd_store(1) = vel_abs(1);
            obj.yd_store(1) = vel_abs(2);
            obj.zx_store(1) = zmp_abs(1);
            obj.zy_store(1) = zmp_abs(2);
            
            obj.QZd = 1;
            obj.QVx = 0;
            obj.QVy = 0;
            obj.QZx = 0;
            obj.QZy = 0;
            
            obj.QCPx = 0;
            obj.QCPy = 10;
              
            obj.footstep = 0;
            obj.first_support_foot = parameters.firstSupportFoot;
            
            obj.vref_x = 0;
            obj.vref_y = 0;
            
            obj.plot_limits = [-1, 1, -1, 1];
            
            obj.footsteps(:,1) = [obj.fxc;obj.fyc];
            obj.orientations(1) = 0;
            
            obj.footsteps_are_preassigned = 0;
            obj.preassigned_footsteps_matrix = [];
            
            % Initialize swing foot
            
            first_left_foot_position = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)]*[0;obj.L/2];
            
            obj.x_l = obj.x + first_left_foot_position(1);
            obj.y_l = obj.y + first_left_foot_position(2);
            obj.xd_l = 0;
            obj.yd_l = 0;
            obj.zx_l = obj.x_l;
            obj.zy_l = obj.y_l;
            
            first_right_foot_position = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)]*[0;-obj.L/2];
            
            obj.x_r = obj.x + first_right_foot_position(1);
            obj.y_r = obj.y + first_right_foot_position(2);
            obj.xd_r = 0;
            obj.yd_r = 0;
            obj.zx_r = obj.x_l;
            obj.zy_r = obj.y_l;
            
            obj.current_support_foot = 'left';
        end
        
        function set_footsteps(obj, preassigned_footsteps_matrix)
%             obj.footsteps_are_preassigned = 1;
            obj.preassigned_footsteps_matrix = preassigned_footsteps_matrix;
        end
        
        function set_vref(obj, vx, vy, ang_vel)
            obj.vref_x = vx;
            obj.vref_y = vy;
            obj.vref_omega = ang_vel;
        end
        
        function set_gains(obj, QZd, QJ, QVx, QVy, QZx, QZy, Qfsx, Qfsy)
            obj.QZd = QZd;
            obj.QJ  = QJ;
            obj.QVx = QVx;
            obj.QVy = QVy;
            obj.QZx = QZx;
            obj.QZy = QZy;
            obj.Qfsx = Qfsx;
            obj.Qfsy = Qfsy;
        end

        function gen_zmp_constraints(obj)
            Ic = zeros(obj.W,obj.W);
            obj.Cc = zeros(obj.W,obj.M);
            
            for i = 1:obj.W
                if not(obj.footstep == 0 && obj.time_discretization_mask(i) <= 2)
                    Ic(i,i) = mod(obj.time_discretization_mask(i),2);
                end
                angle = obj.predicted_rotations(ceil(obj.time_discretization_mask(i)/2));
                rot_cos(i,i) = cos(angle);
                rot_sin(i,i) = sin(angle);
                if obj.time_discretization_mask(i) == 3 || obj.time_discretization_mask(i) == 4
                    obj.Cc(i,1) = 1;
                end
                if obj.time_discretization_mask(i) == 5 || obj.time_discretization_mask(i) == 6
                    obj.Cc(i,2) = 1;
                end
            end
            
            rot = [ rot_cos,-rot_sin; ...
                    rot_sin, rot_cos  ];
            
            AQ_max =  rot'*[   Ic*obj.P , zeros(obj.W, obj.W) , -Ic*obj.Cc , zeros(obj.W, obj.M) ; ...
                               zeros(obj.W, obj.W) , Ic*obj.P , zeros(obj.W, obj.M) , -Ic*obj.Cc ];
                        
            AQ_min = -rot'*[   Ic*obj.P , zeros(obj.W, obj.W) , -Ic*obj.Cc , zeros(obj.W, obj.M) ; ...
                               zeros(obj.W, obj.W) , Ic*obj.P , zeros(obj.W, obj.M) , -Ic*obj.Cc ];
            
            bQ_max = [Ic*obj.p;Ic*obj.p]*obj.w/2 - rot'*[Ic*obj.p*obj.zx;Ic*obj.p*obj.zy] - [obj.pred_zx_foot;obj.pred_zy_foot];
            bQ_min = [Ic*obj.p;Ic*obj.p]*obj.w/2 + rot'*[Ic*obj.p*obj.zx;Ic*obj.p*obj.zy] + [obj.pred_zx_foot;obj.pred_zy_foot];
            
            obj.A_zmp = [AQ_max;AQ_min];
            obj.b_zmp = [bQ_max;bQ_min];
        end

        function gen_footstep_constraints(obj)           
            pfR = zeros(obj.M,1);
            pfL = zeros(obj.M,1);
            if mod(obj.footstep,2) == 0
                for i = 1:obj.M
                    if mod(i,2) == 0
                        pfR(i) = 1;
                    else
                        pfL(i) = 1;
                    end
                end
            else
                for i = 1:obj.M
                    if mod(i,2) == 0
                        pfL(i) = 1;
                    else
                        pfR(i) = 1;
                    end
                end
            end
            
            rot_cos = zeros(obj.M,obj.M);
            rot_sin = zeros(obj.M,obj.M);

            rot_cos(1,1) = cos(0);
            rot_sin(1,1) = sin(0);
            
            for j = 2:obj.M
                % Fill the rotation matrix
                rot_cos(j,j) = cos(obj.predicted_rotations(i));
                rot_sin(j,j) = sin(obj.predicted_rotations(i));
            end

            rot = [ rot_cos,-rot_sin; ...
                    rot_sin, rot_cos  ];
                
            A = eye(obj.M) - [zeros(1,obj.M);eye(obj.M-1),zeros(obj.M-1,1)];
            
            A_fs_max = rot'*[zeros(obj.M,obj.W) , zeros(obj.M,obj.W) , A , zeros(obj.M,obj.M); ...
                             zeros(obj.M,obj.W) , zeros(obj.M,obj.M) , zeros(obj.M,obj.W) , A]; ...
            A_fs_min = rot'*[zeros(obj.M,obj.W) , zeros(obj.M,obj.W) ,-A , zeros(obj.M,obj.M); ...
                             zeros(obj.M,obj.W) , zeros(obj.M,obj.M) , zeros(obj.M,obj.W) ,-A];
            
            obj.A_fs = [A_fs_max;A_fs_min];

            b_fs_max = [ones(obj.M,1)*obj.box_c ; -pfL*obj.box_a + pfR*obj.box_b];
            b_fs_min = [ones(obj.M,1)*obj.box_c ;  pfL*obj.box_b - pfR*obj.box_a];
            
            obj.b_fs = [b_fs_max;b_fs_min];
        end
        
        function gen_stability_constraint_periodic(obj)
            
            obj.Aeq = zeros(2,2*(obj.W+obj.M));
            
            for i = 1:obj.W
                obj.Aeq(1,i) = (1/obj.omega) * (1-exp(-obj.omega*obj.deltas(i))) ...
                               * exp(-obj.omega*sum(obj.deltas(1:i-1))) / (1-exp(-obj.omega*obj.pre_time));
                obj.Aeq(2,obj.W+i) = obj.Aeq(1,i);
            end
            
            %swing_foot_velocity_x = [0 ; obj.pred_zx_foot(2:end) - obj.pred_zx_foot(1:end-1)]/obj.delta;
            %swing_foot_velocity_y = [0 ; obj.pred_zy_foot(2:end) - obj.pred_zy_foot(1:end-1)]/obj.delta;
            swing_foot_velocity_x = zeros(60,1);
            swing_foot_velocity_y = zeros(60,1);
            
            obj.beq = [obj.x + obj.xd/obj.omega - obj.zx; ...
                       obj.y + obj.yd/obj.omega - obj.zy ]...
                       - obj.Aeq * [swing_foot_velocity_x ; swing_foot_velocity_y ; zeros(obj.M,1) ; zeros(obj.M,1)]*0;
        end
        
        function gen_cost_function(obj)
            obj.HQ = obj.H_Zd() + 50000*obj.H_footsteps();
            obj.fQ = obj.f_Zd() + 50000*obj.f_footsteps();
        end
        
        function H = H_Zd(obj)
            HQ1 = obj.QZd*eye(obj.W);
            HQ2 = obj.QZd*eye(obj.W);
            H = blkdiag(HQ1,HQ2,zeros(2*obj.M,2*obj.M));
        end
        
        function f = f_Zd(obj)          
            f = zeros((obj.W+obj.M)*2,1);
        end
        
        function H = H_footsteps(obj)
            
            differenceMatrix = eye(obj.M) - [zeros(1,obj.M);eye(obj.M-1),zeros(obj.M-1,1)];
            
            H = [zeros(obj.W,obj.W),zeros(obj.W,obj.W),zeros(obj.W,obj.M),zeros(obj.W,obj.M); ...
                 zeros(obj.W,obj.W),zeros(obj.W,obj.W),zeros(obj.W,obj.M),zeros(obj.W,obj.M); ...
                 zeros(obj.M,obj.W),zeros(obj.M,obj.W),differenceMatrix,zeros(obj.M,obj.M); ...
                 zeros(obj.M,obj.W),zeros(obj.M,obj.W),zeros(obj.M,obj.M),differenceMatrix];
             
            H = H'*H;
        end
        
        function f = f_footsteps(obj)
            x_integrated = obj.preassigned_footsteps_matrix(1:obj.M,1);
            y_integrated = obj.preassigned_footsteps_matrix(1:obj.M,2);

            differenceMatrix = eye(obj.M) - [zeros(1,obj.M);eye(obj.M-1),zeros(obj.M-1,1)];
            
            f_fs1 = -differenceMatrix'*x_integrated;
            f_fs2 = -differenceMatrix'*y_integrated;
            
            f = [zeros(obj.W,1);zeros(obj.W,1);f_fs1;f_fs2];
        end
        
        function exit_var = cycle(obj,iter)           
            exit_var = false;
            % QP options
            options = optimset('Algorithm','interior-point-convex','Display','off');
            
            % Global iteration and iteration since the start of the last single support phase
            obj.iter = iter;
            obj.iterF = mod(obj.iter,obj.F_samples/2);
            obj.t_curr = obj.iter*obj.delta;
            
            % Construct the time mask
            footstep_mask = [1*ones(1,obj.S_samples),2*ones(1,obj.D_samples)];
            full_mask = [];
            for i = 0:obj.M
                full_mask = [full_mask, footstep_mask + 2*i];
            end
            obj.time_discretization_mask = full_mask(obj.iterF+1:obj.iterF+1+obj.W);
            
            % Construct the integration matrices
            
            obj.p = ones(obj.W,1);
            
            obj.P = [];
            for i = 1:obj.W
                obj.P(i,1:i) = obj.deltas(1:i);
            end
            
            % Position and velocity integration matrices
            
            if obj.iter <= 1
                obj.Pu = [];
                obj.Ps = [];
                obj.Vu = [];
                obj.Vs = [];
                
                for i = 1:obj.W
                    
                    A_power = eye(3);
                    for j = 1:i
                        ch = cosh(obj.omega*obj.deltas(j));
                        sh = sinh(obj.omega*obj.deltas(j));
                        A_upd_j = [ch, sh/obj.omega, 1-ch; obj.omega*sh, ch, -obj.omega*sh; 0, 0, 1];
                        A_power = A_power * A_upd_j;
                    end
                    
                    Ps_newline = A_power(1,:);
                    Vs_newline = A_power(2,:);
                    
                    Pu_newline = zeros(1,obj.W);
                    Vu_newline = zeros(1,obj.W);
                    
                    for j = 1:i
                        A_power = eye(3);
                        for n = j+1:i
                            ch = cosh(obj.omega*obj.deltas(n));
                            sh = sinh(obj.omega*obj.deltas(n));
                            A_upd_n = [ch, sh/obj.omega, 1-ch; obj.omega*sh, ch, -obj.omega*sh; 0, 0, 1];
                            A_power = A_power * A_upd_n;
                            
                        end
                        
                        ch = cosh(obj.omega*obj.deltas(j));
                        sh = sinh(obj.omega*obj.deltas(j));
                        B_upd_j = [obj.deltas(j)-sh/obj.omega; 1-ch; obj.deltas(j)];
                        coeff = A_power*B_upd_j;
                        
                        Pu_newline(1,j) = coeff(1);
                        Vu_newline(1,j) = coeff(2);
                    end
                    
                    obj.Pu = [obj.Pu; Pu_newline];
                    obj.Ps = [obj.Ps; Ps_newline];
                    obj.Vu = [obj.Vu; Vu_newline];
                    obj.Vs = [obj.Vs; Vs_newline];
                end
            end
            
            obj.theta = wrapToPi(obj.theta);
            
            % Current footstep
            obj.footstep = floor(obj.iter/(obj.F_samples/2));

            % Updating time and footsteps
            
            rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
            
            if mod(obj.iter,obj.F_samples/2) == 0 && obj.footstep > 0
                disp(['Changing footstep at ', num2str(obj.t_curr)])
                abs_foot = [obj.fxc;obj.fyc] + rot*[obj.zd(2*obj.W+1);obj.zd(2*obj.W+obj.M+1)];
                obj.fxc = abs_foot(1);
                obj.fyc = abs_foot(2);
                
                obj.footsteps(:,obj.footstep+1) = [obj.fxc;obj.fyc];
                
                obj.theta = obj.theta + obj.predicted_rotations(2); %pi/16;

                obj.orientations(obj.footstep+1) = obj.theta;
                
                new_rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
                
                % Change coordinates relative to new support foot
                pos_rel = new_rot'*rot*([obj.x;obj.y] - [obj.zd(2*obj.W+1);obj.zd(2*obj.W+obj.M+1)]);
                vel_rel = new_rot'*rot*[obj.xd;obj.yd];
                zmp_rel = new_rot'*rot*([obj.zx;obj.zy] - [obj.zd(2*obj.W+1);obj.zd(2*obj.W+obj.M+1)]);

                obj.x = pos_rel(1);
                obj.y = pos_rel(2);
                obj.xd = vel_rel(1);
                obj.yd = vel_rel(2);
                obj.zx = zmp_rel(1);
                obj.zy = zmp_rel(2);
                
                % Update swing foot
                pos_rel = new_rot'*rot*([obj.x_l;obj.y_l] - [obj.zd(2*obj.W+1);obj.zd(2*obj.W+obj.M+1)]);
                vel_rel = new_rot'*rot*[obj.xd_l;obj.yd_l];
                zmp_rel = new_rot'*rot*([obj.zx_l;obj.zy_l] - [obj.zd(2*obj.W+1);obj.zd(2*obj.W+obj.M+1)]);
                
                obj.x_l = pos_rel(1);
                obj.y_l = pos_rel(2);
                obj.xd_l = vel_rel(1);
                obj.yd_l = vel_rel(2);
                obj.zx_l = zmp_rel(1);
                obj.zy_l = zmp_rel(2);
                
                pos_rel = new_rot'*rot*([obj.x_r;obj.y_r] - [obj.zd(2*obj.W+1);obj.zd(2*obj.W+obj.M+1)]);
                vel_rel = new_rot'*rot*[obj.xd_r;obj.yd_r];
                zmp_rel = new_rot'*rot*([obj.zx_r;obj.zy_r] - [obj.zd(2*obj.W+1);obj.zd(2*obj.W+obj.M+1)]);

                obj.x_r = pos_rel(1);
                obj.y_r = pos_rel(2);
                obj.xd_r = vel_rel(1);
                obj.yd_r = vel_rel(2);
                obj.zx_r = zmp_rel(1);
                obj.zy_r = zmp_rel(2);
                
                % If footsteps are preassigned remove past footsteps
                obj.preassigned_footsteps_matrix(1,:) = [];
                
                if strcmp(obj.current_support_foot,'right')
                    obj.current_support_foot = 'left';
                else
                    obj.current_support_foot = 'right';
                end
            end
            
            % Assigning foot rotations
            I_diff_rot = eye(obj.M) - [zeros(1,obj.M);eye(obj.M-1),zeros(obj.M-1,1)];
            p_M = ones(obj.M,1);
            
            angles_integrated = zeros(obj.M,1);
            
            H_rot = (I_diff_rot'*I_diff_rot);
            f_rot = -(I_diff_rot'*angles_integrated);

            optimal_rotations = quadprog(H_rot,f_rot,[],[],[],[],[],[],[],options);
            
            if obj.footstep == 0
                obj.predicted_rotations = [0,0,optimal_rotations(1:obj.M-1)'];
            else
                obj.predicted_rotations = [0,optimal_rotations(1:obj.M)'];
            end
            
            theta_step_max = pi/4;%pi/16;
            
            for i = 2:obj.M+1
                if wrapToPi(angdiff(obj.predicted_rotations(i), obj.predicted_rotations(i-1))) >= theta_step_max
                    obj.predicted_rotations(i) = wrapToPi(obj.predicted_rotations(i-1) + theta_step_max);
                end
                if wrapToPi(angdiff(obj.predicted_rotations(i), obj.predicted_rotations(i-1))) <= -theta_step_max
                    obj.predicted_rotations(i) = wrapToPi(obj.predicted_rotations(i-1) - theta_step_max);
                end
            end
            
            rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
            
            % Swing foot contribution (for now zero)
            obj.pred_zx_foot = 0*ones(obj.W,1);
            obj.pred_zy_foot = zeros(obj.W,1);
            
            % Inequality constraints
            obj.gen_zmp_constraints()
            AQ = obj.A_zmp;
            bQ = obj.b_zmp;
            
            obj.footsteps_are_preassigned = 1;
            % Add bounding boxes on the footsteps if AFP is on
            if obj.footsteps_are_preassigned == 0
                obj.gen_footstep_constraints()
                AQ = [AQ;obj.A_fs];
                bQ = [bQ;obj.b_fs];
            end

            % Stability constraint
            obj.gen_stability_constraint_periodic();
                
%             % Lock the feet in double support
%             if obj.iterF > obj.S_samples
%                 A_lock = zeros(2,2*obj.W+2*obj.M);
%                 A_lock(1,2*obj.W+1) = 1;
%                 A_lock(2,2*obj.W+obj.M+1) = 1;
%                 b_lock = [obj.zd(2*obj.W+1);obj.zd(obj.W*2+obj.M+1)];
%                 obj.Aeq = [obj.Aeq;A_lock];
%                 obj.beq = [obj.beq;b_lock];
%             end
%             
%             % Lock the first footstep in place
%             if obj.footstep == 0 && obj.footsteps_are_preassigned == 0
%                 A_lock = zeros(2,2*obj.W+2*obj.M);
%                 A_lock(1,2*obj.W+1) = 1;
%                 A_lock(2,2*obj.W+obj.M+1) = 1;
%                 b_lock = [0;-obj.L];
%                 obj.Aeq = [obj.Aeq;A_lock];
%                 obj.beq = [obj.beq;b_lock];
%             end
            
            % Cost function
            if size(obj.preassigned_footsteps_matrix,1) == 1
                disp('stop');
                exit_var = true;
                return;
            end
            obj.gen_cost_function()
            
            % Solver
            obj.zd = quadprog(obj.HQ,obj.fQ,AQ,bQ,obj.Aeq,obj.beq,[],[],[],options);
            
            obj.pred_zx_body = obj.p * obj.zx + obj.P*obj.zd(1:obj.W);
            obj.pred_zy_body = obj.p * obj.zy + obj.P*obj.zd(obj.W+1:2*obj.W);
            
            obj.pred_zx = obj.pred_zx_body + obj.pred_zx_foot;
            obj.pred_zy = obj.pred_zy_body + obj.pred_zy_foot;
            
            obj.pred_fs = [obj.zd(2*obj.W+1:2*obj.W+obj.M)';obj.zd(2*obj.W+obj.M+1:2*obj.W+2*obj.M)'];
            
            % State update
            ch = cosh(obj.omega*obj.delta);
            sh = sinh(obj.omega*obj.delta);
            A_upd = [ch, sh/obj.omega, 1-ch; obj.omega*sh, ch, -obj.omega*sh; 0, 0, 1];
            B_upd = [obj.delta-sh/obj.omega; 1-ch; obj.delta];

            x_updated = A_upd*[obj.x; obj.xd; obj.zx] + B_upd*obj.zd(1);
            y_updated = A_upd*[obj.y; obj.yd; obj.zy] + B_upd*obj.zd(obj.W+1);

            obj.x = x_updated(1);
            obj.y = y_updated(1);
            obj.xd = x_updated(2);
            obj.yd = y_updated(2);
            obj.zx = x_updated(3);
            obj.zy = y_updated(3);

            pos_abs = [obj.fxc;obj.fyc] + rot*[obj.x;obj.y];
            vel_abs = rot*[obj.xd;obj.yd];
            zmp_abs = [obj.fxc;obj.fyc] + rot*[obj.zx;obj.zy];
            obj.x_store(iter+1) = pos_abs(1);
            obj.y_store(iter+1) = pos_abs(2);
            obj.xd_store(iter+1) = vel_abs(1);
            obj.yd_store(iter+1) = vel_abs(2);
            obj.zx_store(iter+1) = zmp_abs(1);
            obj.zy_store(iter+1) = zmp_abs(2);
            obj.zxd_store(iter) = obj.zd(1);
            obj.zyd_store(iter) = obj.zd(obj.W+1);
            
            % Update swing foot
            k_sw = 1;
            tf = 0.2;
            t = obj.iterF*obj.delta;
            if t < 0.2
                if mod(obj.footstep,2) == 1
                    obj.xd_l = k_sw*(obj.pred_fs(1,1) - obj.x_l)/(tf-t);
                    obj.yd_l = k_sw*(obj.pred_fs(2,1) - obj.y_l)/(tf-t);
                    obj.x_l = obj.x_l + obj.xd_l*obj.delta;
                    obj.y_l = obj.y_l + obj.yd_l*obj.delta;
                else
                    obj.xd_r = k_sw*(obj.pred_fs(1,1) - obj.x_r)/(tf-t);
                    obj.yd_r = k_sw*(obj.pred_fs(2,1) - obj.y_r)/(tf-t);
                    obj.x_r = obj.x_r + obj.xd_r*obj.delta;
                    obj.y_r = obj.y_r + obj.yd_r*obj.delta;
                end
            end
        end
        
        function set_plot_limits(obj, xmin, xmax, ymin, ymax)
            obj.plot_limits = [xmin, xmax, ymin, ymax];
        end
        
        function set_plot_options(obj, plot_options)
            obj.plot_options = plot_options;
        end
        
        function plot(obj)
            obj.plot_handles = [];
            
            % Rotate the prediction in the world frame
            pred_to_rotate = [obj.pred_zx , obj.pred_zy]';
            rot = [cos(obj.theta),-sin(obj.theta);sin(obj.theta),cos(obj.theta)];
            pred_rotated = rot*pred_to_rotate;
            
            % Plot the predicted ZMP
            if obj.plot_options.plot_pred_zmp
                plot(obj.fxc+pred_rotated(1,:),obj.fyc+pred_rotated(2,:),'m','LineWidth', 2);
            end
            
            % Plot the CoM
            if obj.plot_options.plot_com
                obj.plot_handles(1) = plot(obj.x_store,obj.y_store,'r','LineWidth', 2);
            end
            
            % Plot the ZMP
            if obj.plot_options.plot_zmp
                obj.plot_handles(2) = plot(obj.zx_store,obj.zy_store,'b','LineWidth', 2);
            end
            
            % Plot the footsteps
            if obj.plot_options.plot_footsteps
                for i = 1 : size(obj.footsteps,2)
%                     plot(obj.footsteps(1,i),obj.footsteps(2,i),'ok');
                    obj.plot_constraint(obj.footsteps(1,i), obj.footsteps(2,i), ...
                                        obj.w, obj.w, obj.orientations(i), [0,0,0], 0.0);
                end
            end
            
            % Plot the orientation
            if obj.plot_options.plot_orientation
                quiver(obj.x_store(obj.iter),obj.y_store(obj.iter),0.1*cos(obj.theta),0.1*sin(obj.theta),...
                       'LineWidth',1,'MaxHeadSize',2)
            end

            % Plot predicted footsteps
            rotated_fs = rot*obj.pred_fs;
            rotated_fs = [[0;0],rotated_fs];
            
            for i = 1 : obj.M+1
                % Plot predicted footsteps
                if obj.plot_options.plot_pred_footsteps
                    plot(obj.fxc+rotated_fs(1,i),obj.fyc+rotated_fs(2,i),'ok');
                end
                
                % Plot predicted zmp constraints
                if obj.plot_options.plot_pred_zmp_constraints
                    obj.plot_constraint(obj.fxc+rotated_fs(1,i), obj.fyc+rotated_fs(2,i), ...
                                        obj.w, obj.w, obj.theta+obj.predicted_rotations(i), 'b', 0.1);
                end
                
                % Plot predicted footstep constraints
                if obj.plot_options.plot_pred_footstep_constraints && i< obj.M+1
                    if mod(obj.footstep + i,2) == 0
                        fs_sign = 1;
                    else
                        fs_sign = -1;
                    end
                    obj.plot_constraint(obj.fxc+rotated_fs(1,i) - fs_sign*(obj.box_b/2+obj.box_a/2)*sin(obj.theta+obj.predicted_rotations(i)), ...
                                        obj.fyc+rotated_fs(2,i) + fs_sign*(obj.box_b/2+obj.box_a/2)*cos(obj.theta+obj.predicted_rotations(i)), ...
                                        obj.box_c*2, obj.box_b-obj.box_a, obj.theta+obj.predicted_rotations(i), [204 255 204]/255, 0.7);
                end
            end

            %hold off
            axis equal;
            xlim([obj.plot_limits(1),obj.plot_limits(2)]); ylim([obj.plot_limits(3),obj.plot_limits(4)]);

            %drawnow
        end
        
        function plot_constraint(obj, xc, yc, xw, yw, theta_c, color, alpha)
            center = [xc xc xc xc; yc yc yc yc];
            vertices = [xw/2 xw/2 -xw/2 -xw/2; +yw/2 -yw/2 -yw/2 +yw/2];
            rotation_matrix = [cos(theta_c), -sin(theta_c); sin(theta_c), cos(theta_c)];
            vertices = rotation_matrix*vertices;
            vertices = center + vertices;
            
            point = patch(vertices(1,:),vertices(2,:),color);
            set(point,'FaceAlpha',alpha,'EdgeColor','k','LineWidth',1,'LineStyle','-');
        end
    end
    
end