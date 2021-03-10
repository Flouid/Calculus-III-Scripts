function maxminsaddle(f)

    syms x y;

    % first and second partial derivatives
    fx = diff(f, x);
    fy = diff(f, y);
    fxx = diff(fx, x);
    fxy = diff(fx, y);
    fyy = diff(fy, y);
    % solve fx = 0 and fy = 0 for x and y
    [xcp, ycp] = solve([fx == 0, fy == 0], [x y]);
    % every combination of x and y is a critical point
    crit_points = combvec(xcp', ycp')';
    num_crit_points = length(crit_points(:, 1));

    % EQUATION FOR THE SECOND DERIVATIVE TEST DESCRIMINANT:
    % D(x, y) = fxx(x, y)*fyy(x, y) - fxy(x, y)^2
    % if D > 0 and fxx(x, y) > 0 -> local minimum
    % if D > 0 and fxx(x, y) < 0 -> local minimum
    % if D < 0 -> saddle point
    % if D = 0 -> you're fucked

    % create function for D(x, y)
    D = fxx*fyy - fxy^2;
    % allocate maximum possible side array to mins, maxes, and saddles
    saddle_points_found = zeros([3 num_crit_points]);
    list_minf = zeros([1 num_crit_points]);
    list_maxf = zeros([1 num_crit_points]);
    % create counters to track which index has been assigned
    saddle_points_counter = 1;
    list_minf_counter = 1;
    list_maxf_counter = 1;
    % check the value of D for each critical point
    % calculate D(x, y), fxx(x, y), and f(x, y)
    for p = 1:length(crit_points(:, 1))
        Dp = subs(D, [x y], crit_points(p, :));
        fxxp = subs(fxx, [x y], crit_points(p, :));
        fp = subs(f, [x y], crit_points(p, :));
        % determine what points were found
        if (Dp > 0)
            if (fxxp > 0)
                % a min was found
                if (~all(ismember(fp, list_minf)))
                    % this min has not been seen yet
                    list_minf(list_minf_counter) = fp;
                    list_minf_counter = list_minf_counter + 1;
                    disp('local minimum found:'); disp(fp);
                end
            elseif (fxxp < 0)
                % a max was found
                if (~all(ismember(fp, list_maxf)))
                    % this max has not been seen yet
                    list_maxf(list_maxf_counter) = fp;
                    list_maxf_counter = list_maxf_counter + 1;
                    disp('local maximum found:'); disp(fp);
                end
            end
        elseif (Dp < 0)
            % a saddle point was found
            sp = [crit_points(p, :) fp];
            if (~all(ismember(sp, saddle_points_found)))
                % this saddle point has not been seen yet
                saddle_points_found(:, saddle_points_counter) = sp;
                saddle_points_counter = saddle_points_counter + 1;
                disp('saddle point found:'); disp(sp);
            end
        end
    end
end