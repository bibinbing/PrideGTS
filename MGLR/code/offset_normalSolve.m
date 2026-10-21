function [solution, solver] = offset_normalSolve(normal, rhs, maxConditionNumber, solver)
% Solves a small symmetric positive-definite normal system without inv().

if nargin < 4 || isempty(solver)
    normal = full((normal + normal.') ./ 2);
    if isempty(normal) || size(normal, 1) ~= size(normal, 2) || ...
            any(~isfinite(normal), 'all')
        error('offset:InvalidNormalMatrix', ...
            'The normal matrix must be finite, square, and nonempty.');
    end
    if ~isscalar(maxConditionNumber) || ~isfinite(maxConditionNumber) || ...
            maxConditionNumber <= 1
        error('offset:InvalidMaxConditionNumber', ...
            'maxConditionNumber must be a finite scalar greater than one.');
    end

    reciprocalCondition = rcond(normal);
    [R, flag] = chol(normal);
    if flag ~= 0 || ~isfinite(reciprocalCondition) || ...
            reciprocalCondition < 1 / maxConditionNumber
        error('offset:IllConditionedDesign', ...
            ['The current functional model is rank deficient or exceeds ', ...
            'the configured maximum condition number (estimated condition %.3g).'], ...
            1 ./ max(reciprocalCondition, realmin));
    end
    solver = struct('method', 'chol', 'R', R, ...
        'reciprocalCondition', reciprocalCondition, ...
        'conditionNumber', 1 ./ reciprocalCondition);
end

if size(rhs, 1) ~= size(solver.R, 1)
    error('offset:NormalSolveDimensionMismatch', ...
        'The normal-system right-hand side has an incompatible row count.');
end
solution = solver.R \ (solver.R.' \ rhs);

end
