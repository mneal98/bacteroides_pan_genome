function [h, chi, p] = chi2ind(x, alpha)
    expected = (sum(x,2)*sum(x))/sum(sum(x));
    OE = ((x - expected).^2)./expected;
    chi = sum(sum(OE,2));
    df = (size(x, 1)-1)*(size(x,2)-1);
    p = 1 - chi2cdf(chi,df);
    h = p>=alpha;
end