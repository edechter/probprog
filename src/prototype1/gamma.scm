
(define g  7)
(define p '(0.99999999999980993
            676.5203681218851
            -1259.1392167224028
            771.32342877765313
            -176.61502916214059
            12.507343278686905
            -0.13857109526572012
            9.9843695780195716e-6
            1.5056327351493116e-7))
 
(define (gamma z)
    (if (< (real-part z)  0.5)
        (/ pi ( * (sin (* pi z))  (gamma (- 1 z))))
        ( let ( (z (- z 1))
                (x-init (list-ref p 0)))
          (define (go i sum)
            (if (< i (+ g 2))
                (go (+ i 1) (+ sum (/ (list-ref p i) (+ z i))))
                sum))
          (let ((x (go 1 x-init))
                (t (+ z g 0.5)))
            (* (sqrt (* 2 pi))
               (expt t (+ z 0.5))
               (* x (exp (- t))))))))


