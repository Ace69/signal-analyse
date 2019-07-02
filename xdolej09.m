% Nacte osobni signal
% Vypise vzorkovaci frekvenc, delku ve vzorcich a sekundach  a pocet jednicek a nul
% -------------------------------------------------------------------------------
fprintf('1)\n');
[signal, fs] = audioread('xdolej09.wav'); % nacte signal
fprintf('\t Vzorkovací frekvence Fs = %d kHz\n', fs/1000);
signal=signal'; % prevede signal na radkovy vektor
sample=length(signal);
sec=sample/fs;
bin=sample/16;
fprintf('\t Delka signalu ve vzorcich je %d\n', sample);
fprintf('\t Delka signalu v sekundach je je %f s\n', sec);
fprintf('\t Pocet binarnich symbolu je je %d\n', bin);

% Dekoduje s[n] zpet do binarnich symbolu
%--------------------------------------------------------------------------------
fprintf('2)\n');
j=1;
for i=8:16:bin % dvojity for cyklus pro ulozeni binarnich symbolu
      if signal(i) > 0
          out(j)= 1; 
      else
          out(j)= 0;
       end
       j=j+1;
end

filetext = textread('xdolej09.txt');  % funkce na kontrolu pomoci XOR
%for k=1:1:1999
%    C=xor(filetext(k),out(k));
%    fprintf('\t Xor je: %d\n', C);
%end

time_var = 1/fs; % zajisteni rozmezi, pripraveni na plot
t = 0:time_var:(sample*time_var)-time_var;

plot(t,signal); % plot signalu, upraveni rozsahu
xlim([0 0.02]);
hold on % funkce hold on, aby se puvodni graf neprepsal
stem((1:20)*0.002/2-0.001/2, filetext(1:20)); % vyploceni binarnich symbolu
hold off % funkce hold off, aby se nasledujici vyploceny graf uz prepsal
xlabel('t'); % popis jednotlivych os
ylabel('s[n], symbols');
grid;

% Filtruje a zobrazi nuly a poly
%--------------------------------------------------------------------------------
fprintf('3)\n');

b = [0.0192 -0.0185 -0.0185 0.0192];
a = [1.0000 -2.8870 2.7997 -0.9113];
% tato funkce je okopirovana z referencni funkce ukozmito.m dostupne na strankach predmetu
N = 20; n = 0:N-1; imp = [1 zeros(1,N-1)]; 
h = filter(b,a,imp);  % vyfiltrovani
H = freqz(b,a,256); f=(0:255) / 256 * fs / 2; 
zplane (b,a);       % funkce na plot nul a polu
p = roots(a); 
if (isempty(p) | abs(p) < 1) % overeni stability
  disp('stabilni...')
else
  disp('NESTABILNI !!!')
end
grid;


% Zobrazi modul kmitoctove charakteristiky daneho filtru
%--------------------------------------------------------------------------------
fprintf('4)\n');

freq = freqz(b, a, fs/2);
yGra = abs(freq);
plot(yGra);
%ylim([0 1]);
xlabel('f [Hz]');
grid;

% Zjisteni, o kolik mame posunout
%--------------------------------------------------------------------------------
fprintf('5)\n');
% posunuli jsme o 15 vzorku



% Filtrovani signalu, posunuti, zobrazeni binarnich symbolu a vykresleni puvodniho signalu
%--------------------------------------------------------------------------------
fprintf('6)\n');
ssn = filter(b,a,signal); % filtrovani signalu filtrem danym v zadani

j=1;
for i=8:16:bin % cyklus na dekodovani binarnich symbolu
      if ssn(i) > 0
          outY(j)= 1; 
      else
          outY(j)= 0;
       end
       j=j+1;
end

xssn = circshift(ssn,-15); % posunuti o 15 vzroku
time_var = 1/fs; % zajisteni rozmezi, pripraveni na plot
t = 0:time_var:(sample*time_var)-time_var;
plot(t, signal); % vykresleni puvodniho signalu
hold on
plot(t, ssn); % vykresleni filtrovaneho singalu
plot(t, xssn); % vykresleni fitrlovaneho a posunuteho signalu
stem((1:20)*0.002/2-0.001/2, outY(1:20)); % vykresleni binarnich symbolu
hold off
xlabel('t');
ylabel('s[n], ss[n], ss_s_h_i_f_t_e_d[n], symbols'); % popis
xlim([0 0.02]);
grid;


%--------------------------------------------------------------------------------
fprintf('7)\n');

k=1;
for l=8:16:bin % for cyklus na dekodovani filtrovanyho signalu
      if ssn(l) > 0
          outY(k)= 1; 
      else
          outY(k)= 0;
       end
       k=k+1;
end
pocet=0;
C= xor(out,outY);% xor binarnich symbolu
V=sum(C); % secteni


% Vypociani DFT a vlozeni modulu spekter jak filtrovaneho, tak puvodniho singalu
%--------------------------------------------------------------------------------
fprintf('8)\n');

frequency = (0 : sample / 2 - 1) / sample * fs;
fft_filtered_signal = abs(fft(filter(b, a, signal))); % vypoctiani fast fourier transform z filtrovaneho signalu
fft_filtered_signal = fft_filtered_signal(1 : sample / 2); % pouze do poloviny Fs
subplot(2,1,1); plot(frequency, fft_filtered_signal);
grid;
xlabel('f [Hz]');
legend('Filtered signal');
no_filter= abs(fft(signal)); % vypocitani fft z puvodniho signalu
no_filter = no_filter(1 : sample / 2); % pouze do polovny fs
subplot(2,1,2); plot(frequency, no_filter);
xlabel('f [Hz]');
legend('Non-filtered signal');
grid;


% Odhad funkce hustoty rozdeleni pravdepodobnosti
%--------------------------------------------------------------------------------
fprintf('9)\n');

xHist = hist(signal,50);
plot(xHist);

% Spocitani korelancnich keoficientu od <-50,50>
%--------------------------------------------------------------------------------
fprintf('\n10)\n');
range = (-50 : 50);
R = xcorr(signal,"biased") / sample;
R = R(range + sample);
plot(range , R);
xlabel('range');
grid;


% Spocitani koeff R[0], R[1], R[16]
%--------------------------------------------------------------------------------
fprintf('11)\n');
fprintf('Koeficient R[0] je %f.\n', R(51)); % R[0] lezi uprostred grafu
fprintf('Koeficient R[1] je %f.\n', R(52));
fprintf('Koeficient R[16] je %f.\n', R(68));
 

% Casovy odhad zdruzene funkce hustoty rozdeleni pravdepodobnosti mezi n a n+1 prvekm
% Je zde vyzita funke hist2opt.m, ktere udela vice mene celou ulohu za nas
%--------------------------------------------------------------------------------
fprintf('12)\n');

linSp = linspace(min(signal), max(signal), 100);
[h,p,r] = hist2opt(signal(1:sample-1), signal(2:sample), linSp);
imagesc(-linSp,linSp,p);
axis xy;

% Overeni, jestli integral vyjde 1
% funkce opet prevzata z hist2opt.m
% Funkce hist2opt.m obsahuje funkci check, ktera toto zajisti
%check = sum(sum (p)) * surf; 
%--------------------------------------------------------------------------------
fprintf('13)\n');

% Srovnani hodnoty R[1], vychazi ~0.23, coz je spravne
%Vuzil jsem cast funkce hist2opt.m
% disp(['hist2: check -- 2d integral should be 1 and is ' num2str(check)]); 
%--------------------------------------------------------------------------------
fprintf('14)\n');
disp(r);